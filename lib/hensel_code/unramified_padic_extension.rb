# frozen_string_literal: true

module HenselCode
  # Unramified p-adic extension Q_p(theta), where theta is a root of a monic
  # integer polynomial that is irreducible modulo p.
  #
  # Elements are represented in the fixed power basis 1, theta, ..., theta^(d-1),
  # with d the degree of the defining polynomial. Every element carries its
  # effective p-adic precision: its coefficients are known modulo p^precision.
  # The degree of the basis and the p-adic precision are independent
  # quantities: reduction modulo the defining polynomial and truncation of
  # coefficients modulo p^precision are handled separately and never mixed up.
  #
  # Only the unramified case is supported; a modulus that is not monic or not
  # irreducible modulo p is rejected.
  class UnramifiedPadicExtension
    include Tools

    attr_reader :prime, :modulus, :precision

    def initialize(prime, modulus, precision)
      @prime = prime
      @modulus = modulus
      @precision = precision
      valid_prime?
      valid_precision?
      valid_modulus?
    end

    def ==(other)
      other.is_a?(UnramifiedPadicExtension) &&
        prime == other.prime && modulus == other.modulus && precision == other.precision
    end

    def degree
      modulus.size - 1
    end

    def element(coefficients, precision: nil)
      Element.new(self, coefficients, precision || @precision)
    end

    def zero(precision: nil)
      element(Array.new(degree, 0), precision: precision)
    end

    def one(precision: nil)
      element([1] + Array.new(degree - 1, 0), precision: precision)
    end

    def generator(precision: nil)
      element([0, 1], precision: precision)
    end

    # Reduction modulo the defining polynomial only; coefficients are left
    # untouched (no truncation modulo p^precision).
    def reduce_modulus(coefficients)
      coeffs = coefficients.dup
      while coeffs.size > degree
        coeff = coeffs.pop
        shift = coeffs.size - degree
        modulus[0...degree].each_with_index do |m, i|
          coeffs[shift + i] -= coeff * m
        end
      end
      coeffs
    end

    # Truncation of coefficients modulo p^precision only.
    def truncate(coefficients, precision)
      modulus_p = prime**precision
      coefficients.map { |c| c % modulus_p }
    end

    def reduce_and_truncate(coefficients, precision)
      reduced = reduce_modulus(coefficients)
      reduced += Array.new(degree - reduced.size, 0) if reduced.size < degree
      truncate(reduced, precision)
    end

    def ring_add(coeffs1, coeffs2, precision)
      truncate(coeffs1.zip(coeffs2).map { |a, b| a + b }, precision)
    end

    def ring_sub(coeffs1, coeffs2, precision)
      truncate(coeffs1.zip(coeffs2).map { |a, b| a - b }, precision)
    end

    def ring_mul(coeffs1, coeffs2, precision)
      product = Array.new((2 * degree) - 1, 0)
      coeffs1.each_with_index do |a, i|
        coeffs2.each_with_index do |b, j|
          product[i + j] += a * b
        end
      end
      reduce_and_truncate(product, precision)
    end

    def ring_constant(constant, precision)
      truncate([constant] + Array.new(degree - 1, 0), precision)
    end

    def unit?(coeffs)
      coeffs.any? { |c| (c % prime) != 0 }
    end

    # Inverse of a unit in (Z/p^precision Z)[x]/(modulus) via Hensel/Newton
    # lifting of the inverse over the residue field.
    def ring_inv(coeffs, precision)
      raise NonInvertibleElement, "element is not a unit modulo #{prime}" unless unit?(coeffs)

      inverse = fp_inverse(coeffs)
      inverse += Array.new(degree - inverse.size, 0) if inverse.size < degree
      k = 1
      while k < precision
        k = [2 * k, precision].min
        two = ring_constant(2, k)
        inverse = ring_mul(inverse, ring_sub(two, ring_mul(coeffs, inverse, k), k), k)
      end
      inverse
    end

    # Image of theta under the Frobenius automorphism: the unique root of the
    # modulus congruent to theta^p modulo p, Hensel-lifted to full precision.
    def sigma_theta
      @sigma_theta ||= compute_sigma_theta
    end

    # Frobenius automorphism applied to an element given by its coefficients.
    def frobenius(coeffs, precision)
      result = Array.new(degree, 0)
      power = ring_constant(1, precision)
      coeffs.each do |c|
        result = ring_add(result, truncate(power.map { |pc| pc * c }, precision), precision)
        power = ring_mul(power, sigma_theta, precision)
      end
      result
    end

    def ring_eval(poly, point, precision)
      result = Array.new(degree, 0)
      poly.reverse_each do |coeff|
        result = ring_add(ring_mul(result, point, precision), ring_constant(coeff, precision), precision)
      end
      result
    end

    private

    def valid_prime?
      return if prime.is_a?(Integer) && Prime.prime?(prime)

      raise ArgumentError, "prime must be a prime integer"
    end

    def valid_precision?
      return if precision.is_a?(Integer) && precision.positive?

      raise ArgumentError, "precision must be a positive integer"
    end

    def valid_modulus?
      valid_coefficient_types?
      raise ModulusNotMonic, "modulus must be monic" unless modulus.last == 1
      return if irreducible_mod_p?

      raise ModulusNotUnramified,
            "modulus must be irreducible modulo #{prime}; only unramified extensions are supported"
    end

    def valid_coefficient_types?
      condition = modulus.is_a?(Array) && modulus.size >= 2 && modulus.all? { |c| c.is_a?(Integer) }
      raise ArgumentError, "modulus must be an array of at least two integer coefficients" unless condition
    end

    # Rabin's irreducibility test over the residue field F_p.
    def irreducible_mod_p?
      f = fp_modulus
      x = fp_mod([0, 1], f)
      prime_factors_of(degree).each do |q|
        xpq = x
        (degree / q).times { xpq = fp_powmod(xpq, prime, f) }
        return false unless fp_gcd(fp_sub(xpq, x), f).size <= 1
      end
      xp = x
      degree.times { xp = fp_powmod(xp, prime, f) }
      fp_sub(xp, x).empty?
    end

    def prime_factors_of(number)
      factors = []
      n = number
      candidate = 2
      while candidate * candidate <= n
        if (n % candidate).zero?
          factors << candidate
          n /= candidate while (n % candidate).zero?
        end
        candidate += 1
      end
      factors << n if n > 1
      factors
    end

    def compute_sigma_theta
      s = fp_powmod([0, 1], prime, fp_modulus)
      s += Array.new(degree - s.size, 0) if s.size < degree
      k = 1
      while k < precision
        k = [2 * k, precision].min
        f_of_s = ring_eval(modulus, s, k)
        df_of_s = ring_eval(modulus_derivative, s, k)
        s = ring_sub(s, ring_mul(f_of_s, ring_inv(df_of_s, k), k), k)
      end
      s
    end

    def modulus_derivative
      @modulus_derivative ||= modulus.each_with_index.map { |c, i| c * i }[1..-1]
    end

    # Inverse in the residue field F_p[x]/(modulus) via a^(p^d - 2).
    def fp_inverse(coeffs)
      fp_powmod(fp_trim(coeffs), (prime**degree) - 2, fp_modulus)
    end

    def fp_modulus
      fp_trim(modulus)
    end

    def fp_trim(coeffs)
      trimmed = coeffs.map { |c| c % prime }
      trimmed.pop while trimmed.any? && trimmed.last.zero?
      trimmed
    end

    def fp_sub(coeffs1, coeffs2)
      size = [coeffs1.size, coeffs2.size].max
      fp_trim((0...size).map { |i| (coeffs1[i] || 0) - (coeffs2[i] || 0) })
    end

    def fp_mul(coeffs1, coeffs2)
      return [] if coeffs1.empty? || coeffs2.empty?

      product = Array.new(coeffs1.size + coeffs2.size - 1, 0)
      coeffs1.each_with_index do |a, i|
        coeffs2.each_with_index do |b, j|
          product[i + j] = (product[i + j] + (a * b)) % prime
        end
      end
      fp_trim(product)
    end

    def fp_powmod(base, exponent, mod)
      result = [1]
      b = fp_mod(base, mod)
      e = exponent
      while e.positive?
        result = fp_mod(fp_mul(result, b), mod) if e.odd?
        b = fp_mod(fp_mul(b, b), mod)
        e >>= 1
      end
      result
    end

    def fp_mod(coeffs, divisor)
      fp_divmod(coeffs, divisor).last
    end

    def fp_divmod(coeffs, divisor)
      remainder = fp_trim(coeffs)
      divisor = fp_trim(divisor)
      quotient = Array.new([remainder.size - divisor.size + 1, 0].max, 0)
      inverse_lead = mod_inverse(divisor.last, prime)
      while !remainder.empty? && remainder.size >= divisor.size
        coeff = (remainder.last * inverse_lead) % prime
        shift = remainder.size - divisor.size
        quotient[shift] = coeff
        divisor.each_with_index do |d, i|
          remainder[shift + i] = (remainder[shift + i] - (coeff * d)) % prime
        end
        remainder = fp_trim(remainder)
      end
      [quotient, remainder]
    end

    def fp_gcd(coeffs1, coeffs2)
      a = fp_trim(coeffs1)
      b = fp_trim(coeffs2)
      a, b = b, fp_mod(a, b) until b.empty?
      a
    end

    # Element of the unramified p-adic extension, given by coefficients in the
    # fixed power basis together with its effective p-adic precision.
    class Element
      attr_reader :field, :coefficients, :precision

      def initialize(field, coefficients, precision)
        @field = field
        @precision = precision
        valid_coefficients?(coefficients)
        valid_precision?
        @coefficients = field.reduce_and_truncate(coefficients, precision)
      end

      def +(other)
        compatible?(other)
        precision = effective_precision(other)
        field.element(field.ring_add(coefficients, other.coefficients, precision), precision: precision)
      end

      def -(other)
        compatible?(other)
        precision = effective_precision(other)
        field.element(field.ring_sub(coefficients, other.coefficients, precision), precision: precision)
      end

      def *(other)
        compatible?(other)
        precision = effective_precision(other)
        field.element(field.ring_mul(coefficients, other.coefficients, precision), precision: precision)
      end

      def /(other)
        compatible?(other)
        self * other.inverse
      end

      def inverse
        field.element(field.ring_inv(coefficients, precision), precision: precision)
      end

      def frobenius
        field.element(field.frobenius(coefficients, precision), precision: precision)
      end

      def trace
        result = coefficients
        current = coefficients
        (field.degree - 1).times do
          current = field.frobenius(current, precision)
          result = field.ring_add(result, current, precision)
        end
        field.element(result, precision: precision)
      end

      def norm
        result = coefficients
        current = coefficients
        (field.degree - 1).times do
          current = field.frobenius(current, precision)
          result = field.ring_mul(result, current, precision)
        end
        field.element(result, precision: precision)
      end

      def zero?
        coefficients.all?(&:zero?)
      end

      def unit?
        field.unit?(coefficients)
      end

      def ==(other)
        return false unless other.is_a?(Element) && other.field == field

        precision = effective_precision(other)
        field.truncate(coefficients, precision) == field.truncate(other.coefficients, precision)
      end

      def to_s
        terms = coefficients.map.with_index { |c, i| "#{c}#{theta_power(i)}" }.join(" + ")
        "#{terms} (mod #{field.prime}^#{precision})"
      end

      def inspect
        "<UnramifiedPadicExtension::Element: #{self}>"
      end

      private

      def valid_coefficients?(coefficients)
        condition = coefficients.is_a?(Array) && !coefficients.empty? && coefficients.all? { |c| c.is_a?(Integer) }
        raise ArgumentError, "coefficients must be a non-empty array of integers" unless condition
      end

      def valid_precision?
        condition = precision.is_a?(Integer) && precision.positive? && precision <= field.precision
        raise ArgumentError, "precision must be a positive integer not exceeding the field precision" unless condition
      end

      def compatible?(other)
        return if other.is_a?(Element) && other.field == field

        raise IncompatibleOperandTypes, "operands must be elements of the same unramified p-adic extension"
      end

      def effective_precision(other)
        [precision, other.precision].min
      end

      def theta_power(index)
        ["", "*theta", "*theta^#{index}"][index > 1 ? 2 : index]
      end
    end
  end
end
