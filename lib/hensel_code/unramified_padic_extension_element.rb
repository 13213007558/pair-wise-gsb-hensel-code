# frozen_string_literal: true

module HenselCode
  # Element of an UnramifiedPadicExtension, stored as p^v * (a_0, ..., a_{d-1}).
  class UnramifiedPadicExtensionElement
    attr_reader :field, :valuation, :coefficient_precision, :precision
    protected :field

    def initialize(field, coefficients, valuation, absolute_precision)
      initialize_internal(field, coefficients, valuation,
                          absolute_precision - valuation, false)
    end

    def initialize_internal(field, coefficients, valuation, coefficient_precision, zero)
      @field = field
      @valuation = valuation
      @coefficient_precision = coefficient_precision
      @zero = zero
      set_coefficients(coefficients)
      normalize unless zero
      @precision = @valuation + @coefficient_precision
      @valuation = nil if @zero
    end

    def coefficients
      @coefficients.dup
    end

    alias to_a coefficients

    def zero?
      @zero
    end

    def base?
      @coefficients[1..].all?(&:zero?)
    end

    def base_coefficient
      @coefficients.first
    end

    def +(other)
      same_field!(other)
      return other.dup if zero?
      return dup if other.zero?

      target_valuation = [internal_valuation, other.internal_valuation].min
      common_precision = [precision, other.precision].min - target_valuation
      left = aligned_coefficients(internal_valuation, target_valuation, common_precision)
      right = other.aligned_coefficients(other.internal_valuation, target_valuation, common_precision)
      sum = add_at_precision(left, right, common_precision)
      field.build_element(sum, target_valuation, common_precision)
    end

    def -(other)
      self + (-other)
    end

    def -@
      return dup if zero?

      modulus = prime**@coefficient_precision
      coefficients = @coefficients.map { |coefficient| (-coefficient) % modulus }
      field.build_element(coefficients, internal_valuation, @coefficient_precision)
    end

    def *(other)
      same_field!(other)
      return field.zero if zero? || other.zero?

      common_precision = [@coefficient_precision, other.coefficient_precision].min
      product = field.multiply_unit(at_precision(common_precision),
                                    other.at_precision(common_precision), common_precision)
      field.build_element(product, internal_valuation + other.internal_valuation, common_precision)
    end

    def /(other)
      self * other.inverse
    end

    def inverse
      raise ZeroDivisionError, "cannot invert zero" if zero?

      inverse_coefficients = field.inverse_unit(@coefficients, @coefficient_precision)
      field.build_element(inverse_coefficients, -internal_valuation, @coefficient_precision)
    end

    def frobenius(power = 1)
      return dup if zero?
      unless power.is_a?(Integer) && power >= 0
        raise ArgumentError, "Frobenius power must be a non-negative integer"
      end

      reduced_power = power % field.degree
      result = @coefficients
      reduced_power.times { result = apply_frobenius(result, @coefficient_precision) }
      field.build_element(result, internal_valuation, @coefficient_precision)
    end

    def trace
      return field.zero if zero?

      result = dup
      (1...field.degree).each { |power| result = result + frobenius(power) }
      result
    end

    def norm
      return field.zero if zero?

      result = dup
      (1...field.degree).each { |power| result = result * frobenius(power) }
      result
    end

    def ==(other)
      other.is_a?(self.class) && field.same_configuration?(other.field) &&
        zero? == other.zero? && (zero? || normalized_state == other.normalized_state)
    end

    alias eql? ==

    def hash
      [field.prime, field.precision, field.polynomial, zero?, normalized_state].hash
    end

    def dup
      field.build_element(@coefficients, internal_valuation, @coefficient_precision, zero?)
    end

    def to_s
      return "0" if zero?

      @coefficients.map.with_index do |coefficient, index|
        next if coefficient.zero?
        factor = internal_valuation.zero? ? "" : "#{prime}^#{internal_valuation} * "
        basis = index.zero? ? "1" : "omega^#{index}"
        "#{factor}(#{coefficient}#{basis})"
      end.compact.join(" + ")
    end

    def inspect
      "<UnramifiedPadicExtensionElement: #{self} precision=#{@precision}>"
    end

    protected

    def internal_valuation
      zero? ? nil : @valuation
    end

    def normalized_state
      return [0, field.precision, field.pad_to_degree([0])] if zero?

      [@valuation, @coefficient_precision, @coefficients]
    end

    def at_precision(precision)
      modulus = prime**precision
      @coefficients.map { |coefficient| coefficient % modulus }
    end

    def aligned_coefficients(source_valuation, target_valuation, precision)
      shift = source_valuation - target_valuation
      modulus = prime**precision
      at_precision(precision).map { |coefficient| (coefficient * prime**shift) % modulus }
    end

    def apply_frobenius(coefficients, precision)
      root = field.frobenius_root(precision)
      coefficients.each_with_index.reduce([0]) do |result, (coefficient, index)|
        point = index.zero? ? [1] : repeated_power(root, index, precision)
        term = field.multiply_unit(point, [coefficient], precision)
        field.send(:add, result, term, prime**precision)
      end
    end

    private

    def prime
      field.prime
    end

    def set_coefficients(coefficients)
      validate_coefficients(coefficients)
      unless @coefficient_precision.is_a?(Integer) && @coefficient_precision.positive?
        raise ArgumentError, "precision must be greater than valuation"
      end

      modulus = prime**@coefficient_precision
      padded = field.pad_to_degree(coefficients)
      @coefficients = field.reduce(padded.map { |coefficient| coefficient % modulus }, @coefficient_precision)
    end

    def validate_coefficients(coefficients)
      unless @field.is_a?(UnramifiedPadicExtension)
        raise IncompatibleExtensionFields, "field must be an UnramifiedPadicExtension"
      end
      unless @valuation.is_a?(Integer)
        raise ArgumentError, "valuation must be an integer"
      end
      valid = coefficients.is_a?(Array) && !coefficients.empty? && coefficients.size <= field.degree &&
              coefficients.all? { |coefficient| coefficient.is_a?(Integer) }
      raise ArgumentError, "coefficients must be #{field.degree} integers" unless valid
    end

    def normalize
      while @coefficient_precision > 1 && @coefficients.all? { |coefficient| (coefficient % prime).zero? }
        @coefficient_precision -= 1
        @valuation += 1
        modulus = prime**@coefficient_precision
        @coefficients = @coefficients.map { |coefficient| (coefficient / prime) % modulus }
      end

      mark_zero(@valuation + @coefficient_precision) if @coefficient_precision == 1 &&
                                                        @coefficients.all?(&:zero?)
    end

    def add_at_precision(left, right, precision)
      modulus = prime**precision
      sum = left.map.with_index { |coefficient, index| (coefficient + right[index]) % modulus }
      field.reduce(sum, precision)
    end

    def same_field!(other)
      return if other.is_a?(self.class) && field.same_configuration?(other.field)

      raise IncompatibleExtensionFields, "elements must belong to the same unramified extension"
    end

    def mark_zero(absolute_precision)
      @zero = true
      @valuation = 0
      @coefficient_precision = absolute_precision.positive? ? absolute_precision : field.precision
      @coefficients = field.pad_to_degree([0])
    end

    def repeated_power(base, exponent, precision)
      result = [1]
      power = base
      remaining = exponent
      until remaining.zero?
        result = field.multiply_unit(result, power, precision) if remaining.odd?
        power = field.multiply_unit(power, power, precision)
        remaining >>= 1
      end
      result
    end
  end
end
