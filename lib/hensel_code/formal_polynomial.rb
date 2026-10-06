# frozen_string_literal: true

module HenselCode
  # Formal (dense) polynomial with integer coefficients, stored constant term
  # first. Unlike HenselCode::Polynomial, which represents a p-adic expansion
  # with carry propagation between coefficients, FormalPolynomial is a plain
  # element of Z[x]: arithmetic is exact integer arithmetic, and reductions
  # modulo an integer or computations in F_p[x] are always explicit.
  class FormalPolynomial
    include HenselCode::Tools

    attr_reader :coefficients

    def initialize(coefficients)
      raise ArgumentError, "coefficients can't be nil" if coefficients.nil?
      raise ArgumentError, "coefficients must be an array" unless coefficients.is_a?(Array)
      raise ArgumentError, "coefficients can't be empty" if coefficients.empty?
      raise ArgumentError, "coefficients must be integers" unless coefficients.all? { |c| c.is_a?(Integer) }

      @coefficients = self.class.trim(coefficients)
    end

    def self.trim(coefficients)
      trimmed = coefficients.dup
      trimmed.pop while trimmed.size > 1 && trimmed.last.zero?
      trimmed
    end

    def zero?
      coefficients == [0]
    end

    def degree
      coefficients.size - 1
    end

    def leading_coefficient
      coefficients.last
    end

    def monic?
      !zero? && leading_coefficient == 1
    end

    def ==(other)
      other.is_a?(FormalPolynomial) && coefficients == other.coefficients
    end
    alias eql? ==

    def +(other)
      self.class.new pairwise(other) { |a, b| a + b }
    end

    def -(other)
      self.class.new pairwise(other) { |a, b| a - b }
    end

    def -@
      self.class.new coefficients.map { |c| -c }
    end

    def *(other)
      case other
      when Integer
        self.class.new coefficients.map { |c| c * other }
      when FormalPolynomial
        product = Array.new(degree + other.degree + 1, 0)
        coefficients.each_with_index do |a, i|
          other.coefficients.each_with_index { |b, j| product[i + j] += a * b }
        end
        self.class.new product
      else
        raise WrongHenselCodeInputType, "can't multiply FormalPolynomial by #{other.class}"
      end
    end

    # Coefficient-wise reduction modulo +modulus+.
    def mod(modulus)
      self.class.new coefficients.map { |c| c % modulus }
    end

    # Division with remainder in (Z/modulus Z)[x]. The leading coefficient of
    # +divisor+ must be invertible modulo +modulus+; this always holds over
    # F_p for a nonzero divisor.
    def divmod_mod(divisor, modulus)
      raise ZeroDivisionError, "can't divide by the zero polynomial" if divisor.zero?

      quotient = Array.new([degree - divisor.degree + 1, 1].max, 0)
      remainder = mod(modulus).coefficients
      inverse_leading = mod_inverse(divisor.leading_coefficient % modulus, modulus)
      while remainder != [0] && remainder.size - 1 >= divisor.degree
        shift = remainder.size - 1 - divisor.degree
        coefficient = (remainder.last * inverse_leading) % modulus
        quotient[shift] = coefficient
        divisor.coefficients.each_with_index do |d, i|
          remainder[i + shift] = (remainder[i + shift] - (coefficient * d)) % modulus
        end
        remainder = self.class.trim(remainder)
      end
      [self.class.new(quotient), self.class.new(remainder)]
    end

    # Monic greatest common divisor in F_modulus[x]; +modulus+ must be prime.
    def gcd_mod(other, modulus)
      a = mod(modulus)
      b = other.mod(modulus)
      until b.zero?
        _, r = a.divmod_mod(b, modulus)
        a = b
        b = r
      end
      return a if a.zero? || a.monic?

      (a * mod_inverse(a.leading_coefficient, modulus)).mod(modulus)
    end

    # Inverse of this polynomial modulo (+modulus_polynomial+, +prime+) in
    # F_prime[x]. Raises ZeroDivisionError when the polynomials are not
    # coprime in F_prime[x].
    def inverse_mod(modulus_polynomial, prime)
      r0 = modulus_polynomial.mod(prime)
      r1 = mod(prime)
      t0 = self.class.new([0])
      t1 = self.class.new([1])
      until r1.zero?
        quotient, remainder = r0.divmod_mod(r1, prime)
        r0, r1 = r1, remainder
        t0, t1 = t1, (t0 - (quotient * t1)).mod(prime)
      end
      unless !r0.zero? && r0.degree.zero?
        raise ZeroDivisionError, "#{inspect} has no inverse modulo #{modulus_polynomial.inspect} in F_#{prime}[x]"
      end

      (t0 * mod_inverse(r0.coefficients[0], prime)).mod(prime)
    end

    def to_s
      coefficients.map.with_index do |c, i|
        "#{c}#{polynomial_variable(i)}"
      end.join(" + ")
    end

    def inspect
      "<FormalPolynomial: #{self}>"
    end

    private

    def pairwise(other)
      size = [coefficients.size, other.coefficients.size].max
      a = coefficients + Array.new(size - coefficients.size, 0)
      b = other.coefficients + Array.new(size - other.coefficients.size, 0)
      a.zip(b).map { |x, y| yield x, y }
    end

    def polynomial_variable(index)
      i = index > 1 ? 2 : index
      ["", "x", "x^#{index}"][i]
    end
  end
end
