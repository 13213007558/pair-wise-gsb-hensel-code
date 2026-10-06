# frozen_string_literal: true

module HenselCode
  # Simultaneous Hensel lifting of a whole set of pairwise coprime monic
  # factors of a monic integer polynomial from Z/pZ up to Z/p^N Z.
  #
  # Given a monic f in Z[x] and monic factors g_1, ..., g_r, pairwise
  # coprime modulo p, with f congruent to the product of the g_i modulo p,
  # the factors are lifted together, one exponent at a time, to factors
  # modulo p^N. The factors are only required to be coprime modulo p; they
  # do not need to be irreducible, and f does not need to factor over Z.
  #
  # Every lifting level records the lifted factors (kept in the same order
  # as the input factors, preserving the correspondence) and the residual
  # f - prod(g_i) as an integer FormalPolynomial, which is divisible by the
  # modulus of its level.
  class HenselLift
    # exponent:      k in p^k
    # modulus:       p^exponent
    # factors:       lifted factors modulo p^exponent (same order as input)
    # residual:      f - prod(factors) as an integer FormalPolynomial
    Level = Struct.new(:exponent, :modulus, :factors, :residual)

    attr_reader :polynomial, :prime, :initial_factors, :target_exponent, :start_exponent, :levels

    def self.lift(polynomial, factors, prime, target_exponent, start_exponent: 1)
      new(polynomial, factors, prime, target_exponent, start_exponent: start_exponent).tap(&:lift)
    end

    def initialize(polynomial, factors, prime, target_exponent, start_exponent: 1)
      @polynomial = coerce(polynomial)
      @prime = prime
      @target_exponent = target_exponent
      @start_exponent = start_exponent
      @initial_factors = factors.map { |factor| coerce(factor) }
      @levels = []
      validate!
    end

    def lift
      current = initial_factors.map { |factor| factor.mod(prime**start_exponent) }
      record_level(start_exponent, current)
      (start_exponent...target_exponent).each do |exponent|
        current = lift_step(current, exponent)
        record_level(exponent + 1, current)
      end
      self
    end

    # Lifted factors modulo p^target_exponent, in the same order as the
    # input factors.
    def factors
      levels.last.factors
    end

    # Pairs [input factor (mod p^start_exponent), final lifted factor],
    # preserving the correspondence between each input factor and its lift.
    def correspondence
      levels.first.factors.zip(factors)
    end

    private

    def coerce(object)
      object.is_a?(FormalPolynomial) ? object : FormalPolynomial.new(object)
    end

    def product(factors)
      factors.reduce(FormalPolynomial.one) { |acc, factor| acc * factor }
    end

    def record_level(exponent, current)
      modulus = prime**exponent
      factors = current.map { |factor| factor.mod(modulus) }
      levels << Level.new(exponent, modulus, factors, polynomial - product(factors))
    end

    # One linear lifting step: lifts factors from modulo p^exponent to
    # modulo p^(exponent + 1) by solving, for every i simultaneously,
    #   sum_i delta_i * prod_{j != i} g_j  =  e  (mod p),
    # where e = (f - prod g_j) / p^exponent mod p. Since the g_j are
    # pairwise coprime modulo p, the unique solution with
    # deg(delta_i) < deg(g_i) is delta_i = e * (prod_{j != i} g_j)^-1 mod g_i.
    def lift_step(current, exponent)
      step_modulus = prime**exponent
      residual = polynomial - product(current)
      error = FormalPolynomial.new(residual.coefficients.map { |c| (c / step_modulus) % prime })
      deltas = current.each_index.map do |i|
        others_product = product(current.each_with_index.reject { |_, j| j == i }.map(&:first)).mod(prime)
        inverse = others_product.inverse_mod(current[i].mod(prime), prime)
        _, delta = (error * inverse).mod(prime).divmod_mod(prime, current[i].mod(prime))
        delta
      end
      next_modulus = prime**(exponent + 1)
      current.each_index.map do |i|
        correction = FormalPolynomial.new(deltas[i].coefficients.map { |c| c * step_modulus })
        (current[i] + correction).mod(next_modulus)
      end
    end

    def validate!
      validate_exponents!
      validate_prime!
      raise ArgumentError, "polynomial must be monic" unless polynomial.monic?
      raise ArgumentError, "at least one factor is required" if initial_factors.empty?

      modulus = prime**start_exponent
      initial_factors.each do |factor|
        unless factor.mod(modulus).monic?
          raise ArgumentError, "factors must be monic modulo p^start_exponent"
        end
      end
      validate_degrees!
      validate_product!(modulus)
      validate_coprime!
    end

    def validate_exponents!
      unless start_exponent.is_a?(Integer) && start_exponent >= 1
        raise ArgumentError, "start exponent must be an integer >= 1"
      end
      unless target_exponent.is_a?(Integer) && target_exponent >= start_exponent
        raise ArgumentError, "target exponent must be an integer >= start exponent"
      end
    end

    def validate_prime!
      unless prime.is_a?(Integer) && Prime.prime?(prime)
        raise ArgumentError, "prime must be a prime integer"
      end
    end

    def validate_degrees!
      degree_sum = initial_factors.sum(&:degree)
      unless degree_sum == polynomial.degree
        raise ArgumentError, "factor degrees must sum to the polynomial degree"
      end
    end

    def validate_product!(modulus)
      unless (polynomial - product(initial_factors)).mod(modulus).zero?
        raise ArgumentError, "factors do not multiply to the polynomial modulo p^start_exponent"
      end
    end

    def validate_coprime!
      initial_factors.combination(2) do |factor1, factor2|
        unless factor1.gcd_mod(prime, factor2).degree.zero?
          raise NonCoprimeFactors, "factors #{factor1} and #{factor2} are not coprime modulo #{prime}"
        end
      end
    end
  end
end
