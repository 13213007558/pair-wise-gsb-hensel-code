# frozen_string_literal: true

module HenselCode
  # finite p-adic expansion hensel code class
  class FinitePadicExpansion < PAdicBase
    attr_accessor :polynomial

    def modulus
      prime
    end

    def unit_decomposition
      decomposition = super
      return decomposition if decomposition[:residue].zero?

      decomposition.merge(
        teichmuller_root_digits: p_adic_digits(decomposition[:teichmuller_root], decomposition[:unit_exponent]),
        principal_unit_digits: p_adic_digits(decomposition[:principal_unit], decomposition[:unit_exponent])
      )
    end

    def to_a
      hensel_code
    end

    def to_truncated
      TruncatedFinitePadicExpansion.new(prime, exponent, rational)
    end

    def to_s
      hensel_code.map.with_index do |h, i|
        "#{h}#{polynomial_variable(i)}"
      end.join(" + ")
    end

    def inspect
      "<HenselCode: #{polynomial_form}>"
    end

    def inverse
      new_hensel_code = polynomial.inverse.coefficients
      self.class.new prime, exponent, new_hensel_code
    end

    def replace_hensel_code(new_hensel_code)
      super
      @polynomial = Polynomial.new prime, hensel_code
      self
    end

    private

    def polynomial_form
      to_s
    end

    def evaluate(operation, other)
      new_hensel_code = polynomial.send(operation, other.polynomial).coefficients
      self.class.new prime, exponent, new_hensel_code
    end

    def valid_number?(number)
      if number.is_a?(Rational)
        @rational = number
        @exact_zero_input = number.zero?
      elsif number.is_a?(Array) && number.map(&:class).uniq == [Integer] && number.size == exponent
        @hensel_code = number
        @exact_zero_input = false
        decode
        @polynomial = Polynomial.new prime, hensel_code
      else
        message = "number must be a Rational or an\
        Array of integers of size #{exponent}"
        raise WrongHenselCodeInputType, message
      end
    end

    def valid_hensel_code?(new_hensel_code)
      conditions = [
        new_hensel_code.is_a?(Array),
        new_hensel_code.map(&:class).uniq == [Integer],
        new_hensel_code.size == exponent
      ]
      message = "must be an array of integers of size #{exponent}"
      raise WrongHenselCodeInputType, message unless conditions.uniq == [true]
    end

    def encode
      @hensel_code = rational_to_padic_digits
      @polynomial = Polynomial.new prime, hensel_code
    end

    def decode
      number = 0
      hensel_code.each_with_index { |d, i| number += d * (prime**i) }
      @rational = TruncatedFinitePadicExpansion.new(prime, exponent, number).to_r
    end

    def rational_to_padic_digits
      digits = [rational_to_integer(rational)]
      alpha = rational - digits.last
      (exponent - 1).times do
        alpha = reduce_rational_in_terms_of_prime(alpha)
        digits << rational_to_integer(alpha)
        alpha -= digits.last
      end
      digits
    end

    def reduce_rational_in_terms_of_prime(alpha)
      divisor_numerator = alpha.numerator.gcd(prime)
      divisor_denominator = alpha.denominator.gcd(prime)
      if divisor_numerator != 1
        alpha /= divisor_numerator
      elsif divisor_denominator != 1
        alpha *= divisor_denominator
      end
      alpha
    end

    def rational_to_integer(rat)
      (rat.numerator * mod_inverse(rat.denominator, prime)) % prime
    end

    def unit_decomposition_residue
      hensel_code.each_with_index.sum { |digit, index| digit * (prime**index) }
    end

    def p_adic_digits(number, digit_count)
      digits = []
      digit_count.times do
        digits << (number % prime)
        number /= prime
      end
      digits
    end
  end
end
