# frozen_string_literal: true

module HenselCode
  # Closed p-adic ball represented as center + prime**absolute_precision * Z_p.
  class PAdicBall
    include Tools

    INFINITE_PRECISION = Float::INFINITY

    attr_reader :prime, :center, :valuation, :absolute_precision
    alias precision absolute_precision
    alias exponent absolute_precision

    def initialize(*args, prime: nil, center: nil, valuation: nil,
                   absolute_precision: nil, precision: nil)
      positional_prime, positional_center, positional_valuation,
        positional_precision = args

      @prime = prime || positional_prime
      @center_input = center || positional_center
      @valuation_input = valuation || positional_valuation
      @absolute_precision = absolute_precision || precision || positional_precision

      if args.size == 3
        @valuation_input = nil
        @absolute_precision = positional_valuation
      end

      validate_attributes
      canonicalize
    end

    def self.from_truncated(hensel_code)
      valuation = 0
      residue = hensel_code.to_i
      if residue.zero?
        valuation = hensel_code.exponent
      else
        while (residue % hensel_code.prime).zero?
          residue /= hensel_code.prime
          valuation += 1
        end
      end

      new(
        hensel_code.prime,
        hensel_code.to_i,
        valuation,
        hensel_code.exponent
      )
    end

    def +(other)
      compatible?(other)
      self.class.new(
        prime,
        center + other.center,
        nil,
        [absolute_precision, other.absolute_precision].min
      )
    end

    def -(other)
      compatible?(other)
      self.class.new(
        prime,
        center - other.center,
        nil,
        [absolute_precision, other.absolute_precision].min
      )
    end

    def *(other)
      compatible?(other)
      return self.class.new(prime, 0, INFINITE_PRECISION, INFINITE_PRECISION) if zero? || other.zero?

      result_precision = [
        absolute_precision + other.valuation,
        other.absolute_precision + valuation
      ].min
      self.class.new(prime, center * other.center, nil, result_precision)
    end

    def /(other)
      compatible?(other)
      self * other.inverse
    end

    def inverse
      message = "cannot invert #{inspect}: it contains zero"
      raise ZeroDivisionError, message if contains_zero?

      self.class.new(
        prime,
        Rational(1, 1) / center,
        -valuation,
        absolute_precision - (2 * valuation)
      )
    end

    def -@
      self.class.new(prime, -center, valuation, absolute_precision)
    end

    def +@
      self
    end

    def ==(other)
      other.is_a?(self.class) &&
        prime == other.prime &&
        center == other.center &&
        valuation == other.valuation &&
        absolute_precision == other.absolute_precision
    end
    alias eql? ==

    def hash
      [self.class, prime, center, valuation, absolute_precision].hash
    end

    def cover?(number)
      if number.is_a?(self.class)
        return false if prime != number.prime
        return false if infinite_precision? && number.finite_precision?
        return false if number.absolute_precision < absolute_precision
        return valuation_of(center - number.center) >= absolute_precision
      end

      value = Rational(number)
      return value == center if infinite_precision?

      valuation_of(center - value) >= absolute_precision
    end
    alias include? cover?

    def possible_residues(exponent = absolute_precision + 1)
      message = "integral finite balls require a nonnegative finite precision"
      raise WrongHenselCodeInputType, message unless integral_finite_ball?
      validate_residue_exponent!(exponent)
      residue_modulus = prime**exponent
      count_exponent = [exponent - absolute_precision, 0].max
      count = prime**count_exponent

      Array.new(count) do |index|
        (center + (index * (prime**absolute_precision))) % residue_modulus
      end
    end
    alias residues possible_residues

    def unit_residues(exponent = unit_residue_exponent + 1)
      message = "unit residues require a finite precision"
      raise WrongHenselCodeInputType, message if infinite_precision?
      raise WrongHenselCodeInputType, message unless exponent.is_a?(Integer) && exponent >= unit_residue_exponent

      scaled_center = center / (prime**valuation)
      residue_modulus = prime**exponent
      count = prime**(exponent - unit_residue_exponent)

      Array.new(count) do |index|
        (scaled_center + (index * (prime**unit_residue_exponent))) % residue_modulus
      end
    end

    def zero?
      center.zero? && infinite_precision?
    end
    alias exact_zero? zero?

    def contains_zero?
      zero? || (finite_precision? && valuation >= absolute_precision)
    end

    def uncertain_zero?
      finite_precision? && valuation >= absolute_precision
    end

    def infinite_precision?
      absolute_precision == INFINITE_PRECISION
    end

    def finite_precision?
      !infinite_precision?
    end

    def integral?
      valuation >= 0
    end

    def modulus
      prime**absolute_precision
    end

    def radius
      prime**absolute_precision
    end

    def to_r
      Rational(center)
    end

    def to_s
      return "0" if zero?
      return center.to_s if infinite_precision?

      "#{center} + #{prime}^#{absolute_precision} Z_#{prime}"
    end

    def inspect
      "<HenselCode::PAdicBall: #{self}>"
    end

    private

    def validate_attributes
      validate_prime
      validate_center
      validate_precision
      validate_valuation_input
    end

    def validate_prime
      return if prime.is_a?(Integer) && prime.prime?

      raise WrongHenselCodeInputType, "prime must be a prime Integer"
    end

    def validate_center
      return if @center_input.is_a?(Integer) || @center_input.is_a?(Rational)

      raise WrongHenselCodeInputType, "center must be an Integer or Rational"
    end

    def validate_precision
      return if @absolute_precision == INFINITE_PRECISION
      return if @absolute_precision.is_a?(Integer)

      raise WrongHenselCodeInputType, "absolute precision must be an Integer or infinity"
    end

    def validate_valuation_input
      return if @valuation_input.nil?
      return if @valuation_input == INFINITE_PRECISION
      return if @valuation_input.is_a?(Integer)

      raise WrongHenselCodeInputType, "valuation must be an Integer or infinity"
    end

    def canonicalize
      natural_valuation = valuation_of(@center_input)

      if @absolute_precision == INFINITE_PRECISION
        @valuation = natural_valuation
        @center = normalize_center(Rational(@center_input))
      else
        @valuation = [natural_valuation, @absolute_precision].min
        @center = canonical_center(natural_valuation)
      end

      validate_provided_valuation
    end

    def canonical_center(natural_valuation)
      return normalize_center(Rational(0, 1)) if natural_valuation >= absolute_precision

      residue_exponent = absolute_precision - valuation
      unit = Rational(@center_input) / (prime**valuation)
      residue = modular_rational(unit, prime**residue_exponent)
      normalize_center(Rational(residue) * (prime**valuation))
    end

    def normalize_center(number)
      number.denominator == 1 ? number.numerator : number
    end

    def modular_rational(number, modulus)
      (number.numerator * mod_inverse(number.denominator, modulus)) % modulus
    end

    def validate_provided_valuation
      if zero?
        return if @valuation_input.nil? || @valuation_input == INFINITE_PRECISION
      elsif @valuation_input.nil?
        return
      elsif @valuation_input.is_a?(Integer) && @valuation_input <= valuation
        return
      end

      message = "valuation #{@valuation_input.inspect} is not provable for center #{@center_input} " \
                "at precision #{@absolute_precision}; expected #{valuation.inspect}"
      raise WrongHenselCodeInputType, message
    end

    def compatible?(other)
      message = "#{self} is a #{self.class} while #{other} is a #{other.class}"
      raise IncompatibleOperandTypes, message unless other.is_a?(self.class)

      message = "#{self} has prime #{prime} while #{other} has prime #{other.prime}"
      raise HenselCodesWithDifferentPrimes, message unless prime == other.prime
    end

    def integral_finite_ball?
      !infinite_precision? && integral? && absolute_precision >= 0
    end

    def validate_residue_exponent!(exponent)
      return if exponent.is_a?(Integer) && exponent >= 0

      message = "exponent must be a nonnegative Integer"
      raise WrongHenselCodeInputType, message
    end

    def unit_residue_exponent
      absolute_precision - valuation
    end

    def valuation_of(number)
      return INFINITE_PRECISION if number.zero?

      numerator_valuation = prime_divisor_count(number.numerator)
      denominator_valuation = prime_divisor_count(number.denominator)
      numerator_valuation - denominator_valuation
    end

    def prime_divisor_count(number)
      count = 0
      while (number % prime).zero?
        number /= prime
        count += 1
      end
      count
    end
  end
end
