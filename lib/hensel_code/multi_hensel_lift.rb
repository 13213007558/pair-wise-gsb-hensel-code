# frozen_string_literal: true

module HenselCode
  # Simultaneous (multifactor) Hensel lifting.
  #
  # Given a monic polynomial f in Z[x] and monic polynomials g_1, ..., g_r
  # that are pairwise coprime in F_p[x] with f = g_1 * ... * g_r (mod p),
  # the whole set of factors is lifted synchronously, one p-adic level at a
  # time, to factors G_1, ..., G_r with f = G_1 * ... * G_r (mod p^N).
  #
  # The lifted factors keep the order of the input factors, so G_i is always
  # the lift of g_i, and every G_i stays monic. The result also records, for
  # each level k, the product residual f - G_1 * ... * G_r (an integer
  # polynomial whose coefficients are all divisible by p^k).
  class MultiHenselLift
    # factors:        lifted factors modulo p^exponent, in input order
    # input_factors:  the factors given at initialization, reduced mod p
    # residuals:      hash mapping level k to f - product of factors at level k
    Result = Struct.new(:prime, :exponent, :modulus, :factors, :input_factors, :residuals)

    attr_reader :prime, :polynomial, :factors

    def initialize(prime, polynomial, factors)
      @prime = prime
      @polynomial = polynomial
      @factors = factors
      valid_prime?
      valid_polynomial?
      valid_factors?
      valid_factorization?
      valid_pairwise_coprime?
      precompute_bezout_coefficients
    end

    def lift(exponent)
      unless exponent.is_a?(Integer) && exponent >= 1
        raise ArgumentError, "exponent must be a positive integer"
      end

      lifted = base_factors
      residuals = { 1 => residual(lifted) }
      (1...exponent).each do |k|
        lifted = lift_step(lifted, k)
        residuals[k + 1] = residual(lifted)
      end
      Result.new(prime, exponent, prime**exponent, lifted, factors.map { |g| g.mod(prime) }, residuals)
    end

    private

    attr_reader :base_factors, :bezout_coefficients

    def valid_prime?
      raise ArgumentError, "prime can't be nil" if prime.nil?
      raise ArgumentError, "prime must be an integer" unless prime.is_a?(Integer)
      raise ArgumentError, "prime must be a prime number" unless Prime.prime?(prime)
    end

    def valid_polynomial?
      raise ArgumentError, "polynomial must be a FormalPolynomial" unless polynomial.is_a?(FormalPolynomial)
      raise ArgumentError, "polynomial must be monic" unless polynomial.monic?
    end

    def valid_factors?
      condition = factors.is_a?(Array) && factors.size >= 2 && factors.all? { |g| g.is_a?(FormalPolynomial) }
      raise ArgumentError, "factors must be an array of at least two FormalPolynomial objects" unless condition
      raise ArgumentError, "factors must be monic" unless factors.all?(&:monic?)
    end

    def valid_factorization?
      return if product(factors).mod(prime) == polynomial.mod(prime)

      raise WrongHenselCodeInputType, "the product of the factors must be congruent to the polynomial modulo prime"
    end

    def valid_pairwise_coprime?
      base_factors.combination(2) do |a, b|
        next if a.gcd_mod(b, prime).degree.zero?

        raise NonCoprimePolynomialFactors, "factors #{a.inspect} and #{b.inspect} are not coprime modulo #{prime}"
      end
    end

    # s_i = (product of g_j for j != i)^-1 (mod g_i, p). Because the lifted
    # factors never change modulo p, these Bezout coefficients stay valid at
    # every level of the lifting.
    def precompute_bezout_coefficients
      @bezout_coefficients = base_factors.map.with_index do |gi, i|
        others = base_factors.reject.with_index { |_, j| j == i }
        product(others).inverse_mod(gi, prime)
      end
    end

    def base_factors
      @base_factors ||= factors.map { |g| g.mod(prime) }
    end

    def product(polynomials)
      polynomials.reduce(FormalPolynomial.new([1])) { |acc, g| acc * g }
    end

    def residual(lifted)
      polynomial - product(lifted)
    end

    # Synchronous correction step: from factors correct modulo p^k, compute
    # factors correct modulo p^(k+1). The error polynomial e satisfies
    # f - prod(G_i) = p^k * e, and each factor receives the correction
    # p^k * c_i with c_i = e * s_i (mod g_i, p), so that
    # sum of c_i * prod_{j != i} g_j = e (mod p). Since deg(c_i) < deg(g_i),
    # the leading coefficients are untouched and the factors stay monic.
    def lift_step(lifted, k)
      pk = prime**k
      defect = residual(lifted)
      error = FormalPolynomial.new(defect.coefficients.map { |c| c / pk }).mod(prime)
      lifted.map.with_index do |factor, i|
        correction = (error * bezout_coefficients[i]).divmod_mod(base_factors[i], prime).last
        factor + (correction * pk)
      end
    end
  end
end
