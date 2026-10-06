# frozen_string_literal: true

module HenselCode
  # Formal (carry-free) polynomial with integer coefficients stored
  # little-endian: coefficients[i] is the coefficient of x^i.
  #
  # Unlike HenselCode::Polynomial, which represents a p-adic expansion and
  # propagates carries between consecutive coefficients, FormalPolynomial
  # performs ordinary polynomial arithmetic over the integers, and over Z/mZ
  # through the *_mod helpers. It is therefore the correct representation
  # for elements of F_p[x] and Z[x], e.g. for Hensel lifting of polynomial
  # factorizations.
  class FormalPolynomial
    include HenselCode::Tools

    attr_reader :coefficients

    def initialize(coefficients)
      raise ArgumentError, "coefficients can't be nil" if coefficients.nil?
      raise ArgumentError, "coefficients must be an array" unless coefficients.is_a?(Array)

      integer_coefficients = coefficients.all? { |c| c.is_a?(Integer) }
      raise ArgumentError, "coefficients must be integers" unless integer_coefficients
      raise ArgumentError, "coefficients can't be empty" if coefficients.empty?

      @coefficients = self.class.trim(coefficients)
    end

    def self.trim(coefficients)
      trimmed = coefficients.dup
      trimmed.pop while trimmed.size > 1 && trimmed.last.zero?
      trimmed
    end

    def self.zero
      new([0])
    end

    def self.one
      new([1])
    end

    def degree
      coefficients.size - 1
    end

    def zero?
      coefficients == [0]
    end

    def monic?
      coefficients.last == 1
    end

    def leading_coefficient
      coefficients.last
    end

    def ==(other)
      other.is_a?(FormalPolynomial) && coefficients == other.coefficients
    end

    def +(other)
      self.class.new(zip_coefficients(coefficients, other.coefficients, :+))
    end

    def -(other)
      self.class.new(zip_coefficients(coefficients, other.coefficients, :-))
    end

    def -@
      self.class.new(coefficients.map { |c| -c })
    end

    def *(other)
      product = Array.new(coefficients.size + other.coefficients.size - 1, 0)
      coefficients.each_with_index do |a, i|
        other.coefficients.each_with_index do |b, j|
          product[i + j] += a * b
        end
      end
      self.class.new(product)
    end

    # Exact division over the integers: at every step the leading coefficient
    # of the divisor must divide the leading coefficient of the remainder
    # (always the case when the divisor is monic).
    def divmod(other)
      raise ZeroDivisionError, "division by the zero polynomial" if other.zero?

      quotient = Array.new([degree - other.degree + 1, 1].max, 0)
      remainder = coefficients.dup
      while remainder.size >= other.coefficients.size && remainder != [0]
        shift = remainder.size - other.coefficients.size
        ratio = exact_ratio(remainder.last, other.coefficients.last)
        quotient[shift] = ratio
        other.coefficients.each_with_index { |c, i| remainder[i + shift] -= ratio * c }
        remainder = self.class.trim(remainder)
      end
      [self.class.new(quotient), self.class.new(remainder)]
    end

    # Division over Z/mZ: coefficients are reduced modulo m and the leading
    # coefficient of the divisor must be invertible modulo m.
    def divmod_mod(modulus, other)
      divisor = other.mod(modulus)
      raise ZeroDivisionError, "division by the zero polynomial" if divisor.zero?

      quotient = Array.new([degree - divisor.degree + 1, 1].max, 0)
      remainder = mod(modulus).coefficients
      leading_inverse = mod_inverse(divisor.coefficients.last, modulus)
      while remainder.size >= divisor.coefficients.size && remainder != [0]
        shift = remainder.size - divisor.coefficients.size
        ratio = (remainder.last * leading_inverse) % modulus
        quotient[shift] = ratio
        divisor.coefficients.each_with_index do |c, i|
          remainder[i + shift] = (remainder[i + shift] - (ratio * c)) % modulus
        end
        remainder = self.class.trim(remainder)
      end
      [self.class.new(quotient), self.class.new(remainder)]
    end

    def mod(modulus)
      self.class.new(coefficients.map { |c| c % modulus })
    end

    # Monic gcd of the two polynomials over Z/mZ (m expected to be prime).
    def gcd_mod(modulus, other)
      a = mod(modulus)
      b = other.mod(modulus)
      until b.zero?
        _, remainder = a.divmod_mod(modulus, b)
        a, b = b, remainder
      end
      a.monic_mod(modulus)
    end

    # Extended gcd over Z/mZ: returns [g, s, t] with g monic and
    # s*self + t*other congruent to g modulo m.
    def extended_gcd_mod(modulus, other)
      old_r, r = mod(modulus), other.mod(modulus)
      old_s, s = self.class.one, self.class.zero
      old_t, t = self.class.zero, self.class.one
      until r.zero?
        quotient, remainder = old_r.divmod_mod(modulus, r)
        old_r, r = r, remainder
        old_s, s = s, (old_s - (quotient * s)).mod(modulus)
        old_t, t = t, (old_t - (quotient * t)).mod(modulus)
      end
      inverse = mod_inverse(old_r.coefficients.last, modulus)
      [old_r.monic_mod(modulus), old_s.scale_mod(inverse, modulus), old_t.scale_mod(inverse, modulus)]
    end

    # Inverse of the polynomial modulo another polynomial over Z/mZ.
    def inverse_mod(modulus_polynomial, modulus)
      gcd, s, = extended_gcd_mod(modulus, modulus_polynomial)
      unless gcd.degree.zero?
        raise ArgumentError, "polynomial has no inverse modulo the given polynomial"
      end

      s.mod(modulus)
    end

    def monic_mod(modulus)
      return self if zero?

      scale_mod(mod_inverse(coefficients.last, modulus), modulus)
    end

    def scale_mod(factor, modulus)
      self.class.new(coefficients.map { |c| (c * factor) % modulus })
    end

    def to_s
      terms = coefficients.each_with_index.filter_map do |c, i|
        next if c.zero?

        case i
        when 0 then c.to_s
        when 1 then "#{c}*x"
        else "#{c}*x^#{i}"
        end
      end
      terms.empty? ? "0" : terms.reverse.join(" + ")
    end

    def inspect
      "<FormalPolynomial: #{self}>"
    end

    private

    def zip_coefficients(coefficients1, coefficients2, operation)
      size = [coefficients1.size, coefficients2.size].max
      padded1 = coefficients1 + Array.new(size - coefficients1.size, 0)
      padded2 = coefficients2 + Array.new(size - coefficients2.size, 0)
      padded1.zip(padded2).map { |a, b| a.send(operation, b) }
    end

    def exact_ratio(dividend, divisor)
      raise ZeroDivisionError, "leading coefficient of the divisor is zero" if divisor.zero?

      quotient, remainder = dividend.divmod(divisor)
      unless remainder.zero?
        raise ArgumentError, "inexact division: #{dividend} is not divisible by #{divisor}"
      end

      quotient
    end
  end
end
