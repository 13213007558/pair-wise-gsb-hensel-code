# frozen_string_literal: true

module HenselCode
  # base hensel code class
  class PAdicBase
    include Tools
    include PAdicVerifier

    attr_accessor :prime, :exponent, :rational, :hensel_code, :n
    attr_reader :zero_status
    private :prime=, :exponent=, :rational=, :hensel_code=

    def initialize(prime, exponent, number)
      can_initilize?
      @prime = prime
      @exponent = exponent
      @n = Integer.sqrt(((prime**exponent) - 1) / 2)
      @hensel_code_input = false
      valid_number?(number)
      if @hensel_code_input
        decode
      else
        encode
        decode
      end
      set_zero_status(number)
    end

    def numerator
      rational.numerator
    end

    def denominator
      rational.denominator
    end

    def residue
      padic_integer
    end

    alias to_padic_integer residue

    def valuation
      return nil unless zero_status == :nonzero

      value = padic_integer
      valuation = 0
      while (value % prime).zero?
        value /= prime
        valuation += 1
      end
      valuation
    end

    def nonzero?
      zero_status == :nonzero
    end

    def zero?
      zero_status == :zero
    end

    def zero_mod_prime_power?
      zero_status == :zero_mod_prime_power
    end

    alias known_zero_mod_prime_power? zero_mod_prime_power?

    def unit_decomposition
      k = valuation
      return UnitDecomposition.new(self, zero_status) if k.nil?

      precision = exponent - k
      unit_modulus = prime**precision
      unit = (padic_integer / (prime**k)) % unit_modulus
      root_integer, sign_integer, sign_status = lifted_unit_root(unit, precision)
      principal_integer = principal_unit_residue(unit, root_integer, sign_integer, precision)
      root = build_padic_integer(root_integer, precision)
      principal = build_padic_integer(principal_integer, precision)
      UnitDecomposition.new(self, :nonzero, k, root, principal, sign_integer, sign_status)
    end

    alias multiplicative_unit_decomposition unit_decomposition
    alias unit_factorization unit_decomposition
    alias factorize_unit unit_decomposition
    alias factorize_units unit_decomposition

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
      @hensel_code_input = false
      replace_attribute("rational=", new_rational, 0)
    end

    def replace_hensel_code(new_hensel_code)
      valid_hensel_code?(new_hensel_code)
      @hensel_code_input = true
      replace_attribute("hensel_code=", new_hensel_code, 1)
    end

    private

    def can_initilize?
      message = "#{self.class} can only be inherited."
      raise NonInitializableClass, message if instance_of?(HenselCode::PAdicBase)
    end

    def replace_attribute(attribute, new_value, order)
      send(attribute, new_value)
      if attribute == "hensel_code="
        decode
      else
        order.zero? ? [encode, decode] : [decode, encode]
      end
      set_zero_status(new_value)
      self
    end

    def set_zero_status(number)
      @zero_status = if !padic_integer.zero?
                       :nonzero
                     elsif number.is_a?(Rational) && number.zero?
                       :zero
                     else
                       :zero_mod_prime_power
                     end
    end

    def lifted_unit_root(unit, precision)
      if prime == 2
        if precision < 2
          return [1, nil, :undetermined]
        end

        sign = unit % 4 == 3 ? -1 : 1
        return [sign % (prime**precision), sign, :determined]
      end

      root = unit % prime
      modulus = prime
      (precision - 1).times do
        modulus *= prime
        root = root.pow(prime, modulus)
      end
      [root, 1, :not_applicable]
    end

    def principal_unit_residue(unit, root, sign, precision)
      modulus = prime**precision
      if prime == 2
        sign.nil? ? unit % modulus : (sign * unit) % modulus
      else
        unit * mod_inverse(root, modulus) % modulus
      end
    end

    def padic_integer
      raise NotImplementedError, "#{self.class} must implement #padic_integer"
    end

    def build_padic_integer(_residue, _precision)
      raise NotImplementedError, "#{self.class} must implement #build_padic_integer"
    end
  end
end
