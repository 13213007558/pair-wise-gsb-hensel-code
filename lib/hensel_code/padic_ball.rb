# frozen_string_literal: true

module HenselCode
  # A p-adic ball represented by a center, valuation, and absolute precision.
  class PAdicBall
    include Tools

    INFINITY = Float::INFINITY

    attr_reader :prime, :center, :valuation, :absolute_precision

    alias precision absolute_precision
    alias exponent absolute_precision

    def initialize(prime, center, third_argument = nil, fourth_argument = nil,
                   valuation: nil, absolute_precision: nil, precision: nil, exponent: nil)
      if fourth_argument
        raise WrongHenselCodeInputType, "valuation was supplied twice" if valuation
        unless absolute_precision.nil? && precision.nil? && exponent.nil?
          raise WrongHenselCodeInputType, "absolute precision was supplied twice"
        end
        valuation ||= third_argument
        absolute_precision ||= fourth_argument
      elsif third_argument && (absolute_precision || precision || exponent)
        valuation = third_argument
      elsif third_argument && valuation
        raise WrongHenselCodeInputType, "valuation was supplied twice"
      elsif third_argument
        absolute_precision = third_argument
      end

      if absolute_precision && (precision || exponent)
        raise WrongHenselCodeInputType, "absolute precision was supplied twice"
      end

      absolute_precision ||= precision || exponent
      absolute_precision = INFINITY if absolute_precision.nil? && valuation.nil?
      if valuation && absolute_precision.nil?
        raise WrongHenselCodeInputType, "valuation is required for finite absolute precision"
      end

      @prime = validate_prime(prime)
      @absolute_precision = validate_precision(absolute_precision)
      if !exact? && valuation.is_a?(Integer) && valuation > @absolute_precision
        raise WrongHenselCodeInputType, "valuation cannot exceed absolute precision"
      end
      center_number = validate_center(center)
      @center, inferred_valuation = normalize(center_number, @absolute_precision)

      message = "valuation #{valuation} does not match the provable valuation #{inferred_valuation}"
      supplied_valuation = validate_valuation(valuation) if valuation
      if valuation && !same_valuation?(inferred_valuation, supplied_valuation)
        raise WrongHenselCodeInputType, message
      end

      @valuation = inferred_valuation
    end

    def self.exact(prime, number, valuation: nil)
      new(prime, number, valuation: valuation, absolute_precision: INFINITY)
    end

    def self.zero(prime, absolute_precision)
      new(prime, 0, absolute_precision: absolute_precision)
    end

    def self.from_truncated(hensel_code)
      new(
        hensel_code.prime,
        hensel_code.to_i,
        absolute_precision: hensel_code.exponent
      )
    end

    class << self
      alias from_tfpe from_truncated
    end

    def modulus
      message = "negative absolute precision has no integral modulus"
      raise WrongHenselCodeInputType, message if absolute_precision.negative?
      message = "an exact p-adic value does not have a finite modulus"
      raise WrongHenselCodeInputType, message if exact?

      prime**absolute_precision
    end

    def exact?
      infinite?(absolute_precision)
    end

    def contains_zero?
      center.zero?
    end

    def uncertain?
      contains_zero? && !exact?
    end

    def exact_zero?
      exact? && center.zero?
    end

    alias zero? exact_zero?

    def radius
      message = "exact p-adic values have no nonzero ball radius"
      raise WrongHenselCodeInputType, message if exact?

      prime_power(absolute_precision)
    end

    def +(other)
      valid_operand!(other)
      combine(other, :+)
    end

    def -(other)
      valid_operand!(other)
      combine(other, :-)
    end

    def *(other)
      valid_operand!(other)
      result_precision = [
        absolute_precision + other.valuation,
        valuation + other.absolute_precision
      ].min
      build_ball(center * other.center, result_precision)
    end

    def /(other)
      valid_operand!(other)
      reject_uncertain_inverse!(other)
      return self.class.exact(prime, 0) if exact_zero?
      result_precision = [
        absolute_precision - other.valuation,
        valuation + other.absolute_precision - (2 * other.valuation)
      ].min
      build_ball(center.to_r / other.center, result_precision)
    end

    def inverse
      reject_uncertain_inverse!(self)
      result_precision = exact? ? INFINITY : absolute_precision - (2 * valuation)
      build_ball(Rational(1) / center, result_precision)
    end

    def safe_inverse
      inverse
    rescue ZeroDivisionError
      nil
    end

    def cover?(number)
      value = Rational(number)
      return value == center if exact?

      difference = value - center
      difference.zero? || padic_valuation(difference) >= absolute_precision
    end

    def ==(other)
      other.is_a?(self.class) &&
        prime == other.prime &&
        center == other.center &&
        same_valuation?(valuation, other.valuation) &&
        same_precision?(absolute_precision, other.absolute_precision)
    end

    alias eql? ==

    def hash
      [self.class, prime, center, valuation, absolute_precision].hash
    end

    def to_s
      return center.to_s if exact?

      "#{center} + #{prime}^#{absolute_precision} Z_#{prime}"
    end

    def inspect
      "<PAdicBall: #{self} (v=#{valuation}, N=#{absolute_precision})>"
    end

    private

    def validate_prime(value)
      message = "prime must be a prime Integer"
      raise WrongHenselCodeInputType, message unless value.is_a?(Integer) && Prime.prime?(value)

      value
    end

    def validate_precision(value)
      return INFINITY if infinite?(value)
      return value if value.is_a?(Integer)

      raise WrongHenselCodeInputType, "absolute precision must be an Integer or infinity"
    end

    def validate_valuation(value)
      return INFINITY if infinite?(value)
      return value if value.is_a?(Integer)

      raise WrongHenselCodeInputType, "valuation must be an Integer or infinity"
    end

    def validate_center(value)
      return value if value.is_a?(Integer) || value.is_a?(Rational)

      raise WrongHenselCodeInputType, "center must be an Integer or Rational"
    end

    def validate_operand(other)
      message = "#{other.inspect} is a #{other.class}, not a PAdicBall"
      raise IncompatibleOperandTypes, message unless other.is_a?(PAdicBall)
      return if prime == other.prime

      message = "#{self} has prime #{prime} while #{other} has prime #{other.prime}"
      raise HenselCodesWithDifferentPrimes, message
    end

    alias valid_operand! validate_operand

    def combine(other, operation)
      result_precision = minimum_precision(absolute_precision, other.absolute_precision)
      build_ball(center.send(operation, other.center), result_precision)
    end

    def build_ball(number, result_precision)
      self.class.new(prime, number, absolute_precision: result_precision)
    end

    def reject_uncertain_inverse!(ball)
      return unless ball.contains_zero?

      message = "cannot invert #{ball}: the ball contains zero and its inverse is not determined"
      raise ZeroDivisionError, message
    end

    def normalize(number, precision)
      return [integerize(number), padic_valuation(number)] if infinite?(precision)

      value_valuation = padic_valuation(number)
      return [0, precision] if infinite?(value_valuation) || value_valuation >= precision

      residue_exponent = precision - value_valuation
      unit = number / prime_power(value_valuation)
      modulus = prime**residue_exponent
      unit_residue = (unit.numerator * mod_inverse(unit.denominator, modulus)) % modulus
      [integerize(prime_power(value_valuation) * unit_residue), value_valuation]
    end

    def padic_valuation(number)
      return INFINITY if number.zero?

      valuation = 0
      numerator = number.numerator.abs
      denominator = number.denominator.abs

      while (numerator % prime).zero?
        numerator /= prime
        valuation += 1
      end

      while (denominator % prime).zero?
        denominator /= prime
        valuation -= 1
      end

      valuation
    end

    def prime_power(exponent)
      return Rational(1, prime**(-exponent)) if exponent.negative?

      prime**exponent
    end

    def integerize(number)
      number.denominator == 1 ? number.to_i : number
    end

    def minimum_precision(first, second)
      return second if infinite?(first)
      return first if infinite?(second)

      [first, second].min
    end

    def infinite?(value)
      value.is_a?(Float) && value.infinite? == 1
    end

    def same_valuation?(first, second)
      (infinite?(first) && infinite?(second)) || first == second
    end

    def same_precision?(first, second)
      (infinite?(first) && infinite?(second)) || first == second
    end
  end
end
