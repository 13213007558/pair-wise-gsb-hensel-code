# frozen_string_literal: true

module HenselCode
  # Multivariate Hensel lifting for systems of integer polynomial equations.
  #
  # Given a square system F = (f_1, ..., f_n) of MultivariatePolynomial
  # objects, a prime p, a seed solution modulo p, and a target exponent N,
  # the lifter raises the seed to a solution of F modulo p^N, one p-adic
  # digit at a time:
  #
  #   x_{k+1} = x_k - p^k * J(seed)^{-1} * (F(x_k) / p^k)  (mod p^{k+1})
  #
  # The Jacobian at the seed is verified to be invertible modulo p before
  # lifting. Singular seeds are not expanded into branches; the result
  # clearly reports that lifting is not applicable.
  class MultivariateHensel
    include HenselCode::Tools

    # One lifting step from modulus p^k to p^(k+1).
    Level = Struct.new(:exponent, :modulus, :residuals, :correction, :precision, keyword_init: true)

    # Outcome of a lifting attempt.
    class Result
      attr_reader :solution, :modulus, :precision, :levels, :final_residuals, :reason

      def initialize(applicable:, solution: nil, modulus: nil, precision: 0,
                     levels: [], final_residuals: [], reason: nil)
        @applicable = applicable
        @solution = solution
        @modulus = modulus
        @precision = precision
        @levels = levels
        @final_residuals = final_residuals
        @reason = reason
      end

      def applicable?
        @applicable
      end

      # true when substituting the solution back into every equation
      # gives zero modulo p^precision
      def verified?
        applicable? && final_residuals.all?(&:zero?)
      end
    end

    attr_reader :equations, :prime, :seed, :exponent

    def self.lift(equations:, prime:, seed:, exponent:)
      new(equations: equations, prime: prime, seed: seed, exponent: exponent).lift
    end

    def initialize(equations:, prime:, seed:, exponent:)
      @equations = equations
      @prime = prime
      @seed = seed
      @exponent = exponent
      valid_inputs?
    end

    def lift
      return not_applicable("seed #{seed.inspect} is not a solution modulo #{prime}") unless seed_solution?

      jacobian_inverse = invert_matrix_mod(jacobian_at(seed), prime)
      if jacobian_inverse.nil?
        return not_applicable("Jacobian is singular modulo #{prime}; branch enumeration is not supported")
      end

      current = seed.map { |value| value % prime }
      levels = []
      (1...exponent).each do |k|
        current, level = lift_one_level(current, k, jacobian_inverse)
        levels << level
      end

      final_modulus = prime**exponent
      solution = current.map { |value| value % final_modulus }
      Result.new(
        applicable: true,
        solution: solution,
        modulus: final_modulus,
        precision: exponent,
        levels: levels,
        final_residuals: equations.map { |equation| equation.evaluate_mod(solution, final_modulus) }
      )
    end

    private

    def lift_one_level(point, k, jacobian_inverse)
      modulus = prime**k
      residuals = equations.map { |equation| equation.evaluate(point) }
      rhs = residuals.map { |residual| (-residual / modulus) % prime }
      correction = multiply_matrix_vector_mod(jacobian_inverse, rhs, prime)
      new_point = point.zip(correction).map { |value, digit| value + (digit * modulus) }
      level = Level.new(
        exponent: k,
        modulus: modulus,
        residuals: residuals,
        correction: correction,
        precision: k + 1
      )
      [new_point, level]
    end

    def seed_solution?
      equations.all? { |equation| equation.evaluate_mod(seed, prime).zero? }
    end

    def jacobian_at(point)
      equations.map do |equation|
        (0...arity).map { |index| equation.partial(index).evaluate_mod(point, prime) }
      end
    end

    def arity
      equations.first.arity
    end

    def not_applicable(reason)
      Result.new(applicable: false, reason: reason)
    end

    def invert_matrix_mod(matrix, modulus)
      size = matrix.size
      rows = matrix.each_with_index.map do |row, i|
        row.map { |value| value % modulus } + Array.new(size) { |j| i == j ? 1 : 0 }
      end
      size.times do |column|
        pivot = (column...size).find { |r| (rows[r][column] % modulus).nonzero? }
        return nil if pivot.nil?

        rows[column], rows[pivot] = rows[pivot], rows[column]
        scale = mod_inverse(rows[column][column], modulus)
        rows[column] = rows[column].map { |value| (value * scale) % modulus }
        size.times do |r|
          next if r == column

          factor = rows[r][column]
          next if factor.zero?

          rows[r] = rows[r].each_with_index.map { |value, c| (value - (factor * rows[column][c])) % modulus }
        end
      end
      rows.map { |row| row[size..] }
    end

    def multiply_matrix_vector_mod(matrix, vector, modulus)
      matrix.map do |row|
        row.zip(vector).sum { |a, b| a * b } % modulus
      end
    end

    def valid_inputs?
      valid_equations?
      valid_prime?
      valid_seed?
      valid_exponent?
    end

    def valid_equations?
      raise ArgumentError, "equations can't be nil" if equations.nil?
      unless equations.is_a?(Array) && !equations.empty? && equations.all? { |e| e.is_a?(MultivariatePolynomial) }
        raise ArgumentError, "equations must be a non-empty array of MultivariatePolynomial objects"
      end
      unless equations.all? { |e| e.arity == equations.first.arity }
        raise ArgumentError, "equations must share the same arity"
      end
      return if equations.size == equations.first.arity

      raise ArgumentError, "system must be square: #{equations.size} equations for #{equations.first.arity} variables"
    end

    def valid_prime?
      raise ArgumentError, "prime can't be nil" if prime.nil?
      raise ArgumentError, "prime must be an integer" unless prime.is_a?(Integer)
      raise ArgumentError, "prime must be a prime number" unless Prime.prime?(prime)
    end

    def valid_seed?
      return if seed.is_a?(Array) && seed.size == equations.first.arity && seed.all? { |s| s.is_a?(Integer) }

      raise ArgumentError, "seed must be an array of #{equations.first.arity} integers"
    end

    def valid_exponent?
      raise ArgumentError, "exponent must be an integer" unless exponent.is_a?(Integer)
      raise ArgumentError, "exponent must be greater than or equal to 1" if exponent < 1
    end
  end
end
