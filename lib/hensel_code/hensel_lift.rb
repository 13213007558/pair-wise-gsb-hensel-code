# frozen_string_literal: true

module HenselCode
  # Multivariate Hensel lifting.
  #
  # Given a square system of integer polynomial equations (as structured
  # HenselCode::MultivariatePolynomial objects), a prime p, and a known
  # solution of the system modulo p, lifts the solution level by level to a
  # solution modulo p^N, provided the Jacobian of the system at the seed is
  # invertible modulo p. When the Jacobian is singular, lifting is not
  # applicable and this is reported explicitly instead of attempting to
  # enumerate branches.
  module HenselLift
    extend HenselCode::Tools

    # The outcome of a lifting attempt. When the Jacobian at the seed is
    # singular modulo p, +applicable?+ is false and +status+ is
    # +:singular_jacobian+. Otherwise +solution+ maps each variable to its
    # lifted value modulo +modulus+ (p^N), +levels+ holds one entry per
    # precision level with the residuals and effective precision reached,
    # and +residuals+ holds the exact integer values of every equation at
    # the final solution.
    class LiftResult
      attr_reader :status, :message, :prime, :exponent, :solution, :modulus, :levels, :residuals

      def initialize(status:, message: nil, prime: nil, exponent: nil, solution: nil, modulus: nil, levels: [],
                     residuals: [])
        @status = status
        @message = message
        @prime = prime
        @exponent = exponent
        @solution = solution
        @modulus = modulus
        @levels = levels
        @residuals = residuals
      end

      def applicable?
        status == :ok
      end

      def verified?
        applicable? && residuals.all? { |residual| (residual % modulus).zero? }
      end
    end

    class << self
      def lift(equations:, prime:, exponent:, seed:)
        variables = validate_inputs!(equations, prime, exponent, seed)
        validate_seed!(equations, prime, seed)

        jacobian_inverse = invert_matrix_mod(jacobian_matrix(equations, variables, seed, prime), prime)
        return singular_result(prime, exponent) if jacobian_inverse.nil?

        levels, solution = lift_levels(equations, variables, prime, exponent, seed, jacobian_inverse)
        LiftResult.new(
          status: :ok,
          prime: prime,
          exponent: exponent,
          solution: solution,
          modulus: prime**exponent,
          levels: levels,
          residuals: evaluate_system(equations, solution)
        )
      end

      private

      def singular_result(prime, exponent)
        LiftResult.new(
          status: :singular_jacobian,
          message: "Jacobian is singular modulo #{prime}; Hensel lifting is not applicable to this seed",
          prime: prime,
          exponent: exponent
        )
      end

      def lift_levels(equations, variables, prime, exponent, seed, jacobian_inverse)
        current = variables.map { |name| [name, seed[name] % prime] }.to_h
        levels = [level_entry(equations, current, prime, 1)]
        (1...exponent).each do |k|
          current = lift_one_level(equations, variables, current, jacobian_inverse, prime, k)
          levels << level_entry(equations, current, prime**(k + 1), k + 1)
        end
        [levels, current]
      end

      def level_entry(equations, point, modulus, precision)
        { precision: precision, modulus: modulus, residuals: evaluate_system(equations, point).map { |v| v % modulus } }
      end

      def lift_one_level(equations, variables, current, jacobian_inverse, prime, k)
        factor = prime**k
        rhs = evaluate_system(equations, current).map { |value| ((-value) / factor) % prime }
        correction = jacobian_inverse.map { |row| row.zip(rhs).sum { |a, b| a * b } % prime }
        variables.each_with_index.map do |name, i|
          [name, current[name] + (factor * correction[i])]
        end.to_h
      end

      def evaluate_system(equations, point)
        equations.map { |equation| equation.evaluate(point) }
      end

      def jacobian_matrix(equations, variables, point, prime)
        equations.map do |equation|
          variables.map { |variable| equation.partial(variable).evaluate(point) % prime }
        end
      end

      def invert_matrix_mod(matrix, modulus)
        size = matrix.size
        augmented = matrix.each_with_index.map do |row, i|
          row.map { |v| v % modulus } + Array.new(size) { |j| i == j ? 1 : 0 }
        end
        (0...size).each { |col| return nil unless eliminate_column(augmented, col, size, modulus) }
        augmented.map { |row| row[size..] }
      end

      def eliminate_column(augmented, col, size, modulus)
        pivot = (col...size).find { |row| (augmented[row][col] % modulus) != 0 }
        return false if pivot.nil?

        augmented[col], augmented[pivot] = augmented[pivot], augmented[col]
        inverse = mod_inverse(augmented[col][col] % modulus, modulus)
        augmented[col] = augmented[col].map { |v| (v * inverse) % modulus }
        (0...size).each { |row| eliminate_row(augmented, row, col, modulus) }
        true
      end

      def eliminate_row(augmented, row, col, modulus)
        return if row == col

        factor = augmented[row][col] % modulus
        return if factor.zero?

        augmented[row] = augmented[row].each_with_index.map { |v, j| (v - (factor * augmented[col][j])) % modulus }
      end

      def validate_inputs!(equations, prime, exponent, seed)
        validate_equations!(equations)
        validate_prime!(prime)
        validate_exponent!(exponent)
        variables = equations.first.variables
        validate_seed_shape!(variables, seed)
        variables
      end

      def validate_equations!(equations)
        unless equations.is_a?(Array) && !equations.empty? && equations.all? { |e| e.is_a?(MultivariatePolynomial) }
          raise ArgumentError, "equations must be a non-empty array of HenselCode::MultivariatePolynomial"
        end
        variables = equations.first.variables
        unless equations.all? { |equation| equation.variables == variables }
          raise ArgumentError, "all equations must share the same variables in the same order"
        end
        return if equations.size == variables.size

        raise ArgumentError,
              "system must be square: #{equations.size} equations for #{variables.size} variables"
      end

      def validate_prime!(prime)
        raise ArgumentError, "prime must be an integer" unless prime.is_a?(Integer)
        raise ArgumentError, "#{prime} is not prime" unless Prime.prime?(prime)
      end

      def validate_exponent!(exponent)
        return if exponent.is_a?(Integer) && exponent >= 1

        raise ArgumentError, "exponent must be an integer greater than or equal to 1"
      end

      def validate_seed_shape!(variables, seed)
        raise ArgumentError, "seed must be a hash" unless seed.is_a?(Hash)

        missing = variables.reject { |name| seed[name].is_a?(Integer) }
        return if missing.empty?

        raise ArgumentError, "seed must assign an integer to each variable; missing: #{missing.inspect}"
      end

      def validate_seed!(equations, prime, seed)
        residuals = evaluate_system(equations, seed).map { |value| value % prime }
        return if residuals.all?(&:zero?)

        raise ArgumentError, "seed is not a solution of the system modulo #{prime}"
      end
    end
  end
end
