# frozen_string_literal: true

module HenselCode
  # base hensel code class
  class PAdicBase
    include Tools
    include PAdicVerifier

    attr_accessor :prime, :exponent, :rational, :hensel_code, :n
    private :prime=, :exponent=, :rational=, :hensel_code=

    def initialize(prime, exponent, number)
      can_initilize?
      @prime = prime
      @exponent = exponent
      @n = Integer.sqrt(((prime**exponent) - 1) / 2)
      @exact_zero_input = false
      valid_number?(number)
      if number.is_a?(Rational)
        encode
        decode
      end
    end

    def zero_status
      unit_decomposition[:status]
    end

    def exact_zero?
      unit_decomposition[:exact_zero]
    end

    def zero_modulo_prime_power?
      unit_decomposition[:zero_modulo_prime_power]
    end

    def unit_decomposition
      PAdicUnitDecomposition.decompose(
        prime,
        exponent,
        unit_decomposition_residue,
        @exact_zero_input
      )
    end

    def numerator
      rational.numerator
    end

    def denominator
      rational.denominator
    end

    def to_r
      decode
      rational
    end

    def +(other)
      valid?(other)
      evaluate("+", other)
    end

    def -(other)
      valid?(other)
      evaluate("-", other)
    end

    def *(other)
      valid?(other)
      evaluate("*", other)
    end

    def /(other)
      valid?(other)
      evaluate("/", other)
    end

    def replace_prime(new_prime)
      replace_attribute("prime=", new_prime, 0)
    end

    def replace_exponent(new_exponent)
      replace_attribute("exponent=", new_exponent, 0)
    end

    def replace_rational(new_rational)
      @exact_zero_input = new_rational.zero?
      replace_attribute("rational=", new_rational, 0)
    end

    def replace_hensel_code(new_hensel_code)
      valid_hensel_code?(new_hensel_code)
      @exact_zero_input = false
      send("hensel_code=", new_hensel_code)
      decode
      self
    end

    private

    def can_initilize?
      message = "#{self.class} can only be inherited."
      raise NonInitializableClass, message if instance_of?(HenselCode::PAdicBase)
    end

    def unit_decomposition_residue
      hensel_code % (prime**exponent)
    end

    def replace_attribute(attribute, new_value, order)
      send(attribute, new_value)
      order.zero? ? [encode, decode] : [decode, encode]
      self
    end
  end
end
