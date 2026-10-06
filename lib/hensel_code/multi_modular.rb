# frozen_string_literal: true

module HenselCode
  # Multi-modular exact computation of matrix determinants.
  #
  # The determinant is evaluated independently modulo several primes and the
  # signed integer (or rational) determinant is recovered with the Chinese
  # Remainder Theorem once the product of the moduli exceeds twice an explicit
  # Hadamard bound, which guarantees a unique representation in the symmetric
  # residue system. A prime dividing the determinant is a normal event: the
  # residue 0 modulo that prime is simply combined through CRT and never
  # treated as singularity of the integer matrix.
  module MultiModular
    include Tools

    # Result of a multi-modular determinant computation, carrying the exact
    # determinant together with the data that certifies its uniqueness.
    Result = Struct.new(:determinant, :bound, :modulus, :primes, :remainders, keyword_init: true)

    # Explicit Hadamard bound for an integer matrix: |det| <= prod of the
    # Euclidean row norms, computed exactly as floor(sqrt(prod of squared
    # row norms)).
    def hadamard_bound(matrix)
      validate_matrix!(matrix)
      sum_of_squares = matrix.map { |row| row.sum { |entry| entry * entry } }
      Integer.sqrt(sum_of_squares.inject(:*))
    end

    # Determinant of an integer matrix modulo a prime via Gaussian elimination
    # with partial pivoting over the prime field.
    def determinant_mod(matrix, prime)
      validate_matrix!(matrix)
      size = matrix.size
      rows = matrix.map { |row| row.map { |entry| entry % prime } }
      determinant = 1
      size.times do |col|
        pivot = (col...size).find { |i| rows[i][col] != 0 }
        return 0 if pivot.nil?

        unless pivot == col
          rows[col], rows[pivot] = rows[pivot], rows[col]
          determinant = -determinant
        end
        determinant = (determinant * rows[col][col]) % prime
        eliminate_column(rows, col, mod_inverse(rows[col][col], prime), prime)
      end
      determinant
    end

    # Exact determinant of a square matrix of Integers and/or Rationals.
    #
    # Rational entries are cleared row-wise by the exact lcm of the row
    # denominators and the final result is rescaled by the exact product of
    # those lcms. Primes are taken from the +primes+ enumerable when given,
    # otherwise from all primes >= +min_prime+. CRT accumulation stops as soon
    # as the modulus product exceeds twice the Hadamard bound.
    def multi_modular_determinant(matrix, primes: nil, min_prime: 2**16)
      integer_matrix, scale = clear_denominators(matrix)
      bound = hadamard_bound(integer_matrix)
      used_primes = []
      remainders = []
      modulus = 1
      (primes || Prime.each).each do |prime|
        next if primes.nil? && prime < min_prime

        used_primes << prime
        remainders << determinant_mod(integer_matrix, prime)
        modulus *= prime
        break if modulus > 2 * bound
      end
      combined = crt(used_primes, remainders)
      signed = combined > modulus / 2 ? combined - modulus : combined
      determinant = scale == 1 ? signed : Rational(signed, scale)
      Result.new(determinant: determinant, bound: bound, modulus: modulus,
                 primes: used_primes, remainders: remainders)
    end

    private

    def validate_matrix!(matrix)
      raise NonSquareMatrix, "the matrix must not be empty" if matrix.empty?

      size = matrix.size
      unless matrix.all? { |row| row.size == size }
        raise NonSquareMatrix, "the matrix must be square (got #{size} rows with different column counts)"
      end
    end

    def eliminate_column(rows, col, pivot_inverse, prime)
      size = rows.size
      ((col + 1)...size).each do |i|
        factor = (rows[i][col] * pivot_inverse) % prime
        next if factor.zero?

        (col...size).each { |j| rows[i][j] = (rows[i][j] - (factor * rows[col][j])) % prime }
      end
    end

    def clear_denominators(matrix)
      validate_matrix!(matrix)
      row_scales = matrix.map { |row| row.map(&:denominator).inject(:lcm) }
      scaled = matrix.each_with_index.map do |row, i|
        row.map { |entry| (entry * row_scales[i]).to_i }
      end
      [scaled, row_scales.inject(:*)]
    end
  end
end
