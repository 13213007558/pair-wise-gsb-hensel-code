# frozen_string_literal: true

module HenselCode
  # Formal polynomial in a variable independent from the p variable used by Polynomial.
  class FormalPolynomial
    attr_reader :prime, :precision, :coefficients

    def initialize(prime, precision, coefficients)
      validate_base!(prime, precision, coefficients)
      @prime = prime
      @precision = precision
      @coefficients = coefficients.map { |coefficient| normalize(coefficient) }
    end

    def self.from_padic_coefficients(prime, precision, coefficients)
      unless precision.is_a?(Integer) && precision.positive?
        raise ArgumentError, "precision must be positive"
      end

      integer_coefficients = coefficients.map do |coefficient|
        if coefficient.is_a?(Array)
          valid_digits = coefficient.size == precision &&
                         coefficient.all? { |digit| digit.is_a?(Integer) && digit >= 0 && digit < prime }
          unless valid_digits
            raise ArgumentError, "p-adic digit arrays must contain #{precision} integers"
          end

          coefficient.each_with_index.sum { |digit, index| digit * (prime**index) }
        elsif coefficient.respond_to?(:hensel_code) && coefficient.hensel_code.is_a?(Array)
          unless coefficient.prime == prime && coefficient.exponent == precision
            raise ArgumentError, "Hensel code must use the same prime and precision"
          end

          coefficient.hensel_code.each_with_index.sum { |digit, index| digit * (prime**index) }
        elsif coefficient.respond_to?(:hensel_code)
          unless coefficient.prime == prime && coefficient.exponent == precision
            raise ArgumentError, "Hensel code must use the same prime and precision"
          end

          coefficient.hensel_code
        else
          coefficient
        end
      end
      new(prime, precision, integer_coefficients)
    end

    def degree_bound
      coefficients.size - 1
    end

    def degree
      coefficients.rindex { |coefficient| !coefficient.zero? }
    end

    alias known_degree degree

    def leading_residue
      coefficients.last
    end

    def [](index)
      return 0 if index.negative? || index >= coefficients.size

      coefficients[index]
    end

    def +(other)
      compatible?(other)
      self.class.new(prime, precision, combine(other, 1))
    end

    def -(other)
      compatible?(other)
      self.class.new(prime, precision, combine(other, -1))
    end

    def *(other)
      compatible?(other)
      return self.class.new(prime, precision, []) if coefficients.empty? || other.coefficients.empty?

      products = Array.new(coefficients.size + other.coefficients.size - 1, 0)
      coefficients.each_with_index do |left, i|
        other.coefficients.each_with_index do |right, j|
          products[i + j] += left * right
        end
      end
      self.class.new(prime, precision, products)
    end

    def coefficients_with_zeros(size = coefficients.size)
      values = Array.new(size, 0)
      coefficients.each_with_index { |value, index| values[index] = value }
      values
    end

    def valuation_report
      coefficients.map.with_index do |coefficient, index|
        residue = PadicResidue.new(prime, precision, coefficient)
        {
          degree: index,
          residue: residue.residue,
          valuation: residue.valuation_determined? ? residue.valuation : nil,
          valuation_lower_bound: residue.valuation,
          valuation_determined: residue.valuation_determined?,
          known_unit_precision: residue.known_precision,
          unresolved_precision: residue.unresolved_precision
        }
      end
    end

    def to_h
      {
        prime: prime,
        precision: precision,
        degree_bound: degree_bound,
        degree: degree,
        coefficients: coefficients.dup
      }
    end

    def to_s(variable = "x")
      return "0" if coefficients.empty?

      coefficients.map.with_index do |coefficient, degree|
        next if coefficient.zero?

        case degree
        when 0 then coefficient.to_s
        when 1 then "#{coefficient}#{variable}"
        else "#{coefficient}#{variable}^#{degree}"
        end
      end.compact.join(" + ")
    end

    def inspect
      "<FormalPolynomial: #{self}>"
    end

    private

    def validate_base!(prime, precision, coefficients)
      raise ArgumentError, "prime must be an integer" unless prime.is_a?(Integer)
      raise ArgumentError, "prime must be positive" unless prime.positive?
      raise ArgumentError, "precision must be positive" unless precision.positive?
      unless coefficients.is_a?(Array) && coefficients.all? { |value| value.is_a?(Integer) }
        raise ArgumentError, "coefficients must be an array of integers"
      end
    end

    def normalize(coefficient)
      coefficient % (prime**precision)
    end

    def compatible?(other)
      unless other.is_a?(self.class) && other.prime == prime && other.precision == precision
        raise ArgumentError, "polynomials must have the same prime and precision"
      end
    end

    def combine(other, sign)
      size = [coefficients.size, other.coefficients.size].max
      coefficients_with_zeros(size).map.with_index do |value, index|
        value + (sign * other[index])
      end
    end
  end
end
