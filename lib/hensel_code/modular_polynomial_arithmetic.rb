# frozen_string_literal: true

module HenselCode
  # Arithmetic for polynomials over Z/mZ and finite fields F_p.
  module ModularPolynomialArithmetic
    module_function

    def trim(polynomial)
      result = polynomial.dup
      result.pop while result.size > 1 && result.last.zero?
      result
    end

    def reduce_coefficients(polynomial, modulus)
      polynomial.map { |coefficient| coefficient % modulus }
    end

    def add(left, right, modulus)
      length = [left.size, right.size].max
      padded_left = left.dup.fill(0, left.size...length)
      padded_right = right.dup.fill(0, right.size...length)
      padded_left.zip(padded_right).map { |a, b| (a + b) % modulus }
    end

    def negate(polynomial, modulus)
      polynomial.map { |coefficient| (-coefficient) % modulus }
    end

    def subtract(left, right, modulus)
      add(left, negate(right, modulus), modulus)
    end

    def multiply(left, right, modulus)
      product = Array.new(left.size + right.size - 1, 0)
      left.each_with_index do |left_coefficient, i|
        right.each_with_index do |right_coefficient, j|
          product[i + j] += left_coefficient * right_coefficient
        end
      end
      product.map { |coefficient| coefficient % modulus }
    end

    def divmod(numerator, denominator, modulus)
      divisor = trim(reduce_coefficients(denominator, modulus))
      raise ZeroDivisionError, "polynomial divisor is zero" if divisor == [0]

      remainder = trim(reduce_coefficients(numerator, modulus))
      quotient = []
      inverse_leading_coefficient = inverse_integer(divisor.last, modulus)

      while remainder.size >= divisor.size && remainder != [0]
        shift = remainder.size - divisor.size
        factor = remainder.last * inverse_leading_coefficient % modulus
        quotient[shift] = factor.to_i
        divisor.each_with_index do |coefficient, i|
          remainder[shift + i] -= factor * coefficient
          remainder[shift + i] %= modulus
        end
        remainder = trim(remainder)
      end

      normalized_quotient = quotient.map { |coefficient| coefficient.to_i % modulus }
      [trim(normalized_quotient.empty? ? [0] : normalized_quotient), remainder]
    end

    def modulo(polynomial, modulus_polynomial, modulus)
      divmod(polynomial, modulus_polynomial, modulus).last
    end

    def multiply_mod(left, right, modulus_polynomial, modulus)
      product = multiply(left, right, modulus)
      modulo(product, modulus_polynomial, modulus)
    end

    def power_mod(base, exponent, modulus_polynomial, modulus)
      result = [1]
      power = reduce_coefficients(base, modulus)

      while exponent.positive?
        result = multiply_mod(result, power, modulus_polynomial, modulus) if exponent.odd?
        power = multiply_mod(power, power, modulus_polynomial, modulus)
        exponent >>= 1
      end
      result
    end

    def gcd(left, right, modulus_polynomial, modulus)
      left, right = right, modulo(left, right, modulus) until right == [0]
      normalize_monic(left, modulus)
    end

    def bezout(left, right, modulus)
      old_r = trim(reduce_coefficients(left, modulus))
      r = trim(reduce_coefficients(right, modulus))
      old_s = [1]
      polynomial_s = [0]
      old_t = [0]
      polynomial_t = [1]

      until r == [0]
        quotient, remainder = divmod(old_r, r, modulus)
        old_r, r = r, remainder
        old_s, polynomial_s = polynomial_s, subtract(old_s, multiply(quotient, polynomial_s, modulus), modulus)
        old_t, polynomial_t = polynomial_t, subtract(old_t, multiply(quotient, polynomial_t, modulus), modulus)
      end

      inverse = inverse_integer(old_r.last, modulus)
      [normalize_monic(old_r, modulus), old_s.map { |c| c * inverse % modulus },
       old_t.map { |c| c * inverse % modulus }]
    end

    def inverse_mod(polynomial, modulus_polynomial, modulus)
      _, left, _ = bezout(modulo(polynomial, modulus_polynomial, modulus), modulus_polynomial, modulus)
      inverse = modulo(left, modulus_polynomial, modulus)
      unless invertible?(inverse, polynomial, modulus_polynomial, modulus)
        raise ZeroDivisionError, "polynomial is not invertible modulo #{modulus_polynomial}"
      end

      inverse
    end

    def irreducible_over_prime_field?(polynomial, prime)
      f = trim(reduce_coefficients(polynomial, prime))
      degree = f.size - 1
      return false if degree <= 0 || f.last != 1
      return true if degree == 1
      return false if f.first.zero?

      x = [0, 1]
      Prime.prime_division(degree).map(&:first).uniq.each do |prime_factor|
        target_exponent = degree / prime_factor
        power = frobenius_power(x, 0, target_exponent, f, prime)
        next if gcd(subtract(power, x, prime), f, f, prime) == [1]

        return false
      end

      power = frobenius_power(x, 0, degree, f, prime)
      trim(subtract(power, x, prime)) == [0]
    end

    def normalize_monic(polynomial, modulus)
      trimmed = trim(polynomial)
      return trimmed if trimmed == [0]

      inverse = inverse_integer(trimmed.last, modulus)
      trimmed.map { |coefficient| coefficient * inverse % modulus }
    end

    def inverse_integer(number, modulus)
      number = number % modulus
      remainder = modulus
      coefficient = 0
      previous_remainder = number
      previous_coefficient = 1

      until remainder.zero?
        quotient = previous_remainder / remainder
        previous_remainder, remainder = remainder, previous_remainder - (quotient * remainder)
        previous_coefficient, coefficient = coefficient, previous_coefficient - (quotient * coefficient)
      end

      gcd = previous_remainder
      raise ZeroDivisionError, "#{number} has no inverse modulo #{modulus}" unless gcd == 1

      previous_coefficient % modulus
    end

    def frobenius_power(power, exponent, target_exponent, modulus_polynomial, prime)
      (target_exponent - exponent).times do
        power = power_mod(power, prime, modulus_polynomial, prime)
      end
      power
    end

    def invertible?(inverse, polynomial, modulus_polynomial, modulus)
      multiply_mod(inverse, polynomial, modulus_polynomial, modulus) == [1]
    end
  end
end
