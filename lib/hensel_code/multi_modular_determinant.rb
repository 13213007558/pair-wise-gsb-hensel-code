# frozen_string_literal: true

module HenselCode
  # Exact determinant of integer and rational matrices via multi-modular
  # evaluation.
  #
  # The determinant is evaluated independently modulo each prime from a
  # sequence of primes, and the residues are combined with the Chinese
  # remainder theorem. An explicit Hadamard bound on the absolute value of
  # the determinant decides when the product of the primes is large enough
  # to uniquely recover the signed integer determinant.
  #
  # A prime dividing the determinant is a normal event: the residue modulo
  # that prime is simply zero and takes part in the CRT combination like
  # any other residue, so matrices with a zero determinant (globally, or
  # only modulo the first few primes) are handled correctly.
  #
  # Rational entries are supported by clearing denominators row-wise with
  # exact integer arithmetic; the final result is rescaled by the exact
  # product of the cleared denominators.
  class MultiModularDeterminant
    include HenselCode::Tools

    # @param matrix [Array<Array<Integer, Rational>>] square matrix with
    #   Integer or Rational entries
    # @param primes [Array<Integer>, Enumerator, nil] primes to use for the
    #   modular evaluations, in order; defaults to all primes in ascending
    #   order
    def initialize(matrix, primes: nil)
      @matrix = validate_matrix(matrix)
      @prime_source = primes.nil? ? Prime.each : primes.each
      @primes = []
      @residues = []
      @modulus = 1
      @computed = false
    end

    def self.compute(matrix, primes: nil)
      new(matrix, primes: primes).determinant
    end

    # Exact determinant of the input matrix. Returns an Integer for integer
    # matrices and a Rational for matrices with rational entries.
    def determinant
      compute unless @computed
      @determinant
    end

    def bound
      compute unless @computed
      @bound
    end

    def primes
      compute unless @computed
      @primes
    end

    def residues
      compute unless @computed
      @residues
    end

    def modulus
      compute unless @computed
      @modulus
    end

    def scale
      compute unless @computed
      @scale
    end

    private

    def compute
      integer_matrix, @scale = clear_denominators(@matrix)
      @bound = hadamard_bound(integer_matrix)
      @determinant =
        if @bound.zero?
          0
        else
          recover_determinant(integer_matrix)
        end
      @determinant = Rational(@determinant, @scale) unless @scale == 1
      @computed = true
    end

    def recover_determinant(integer_matrix)
      loop do
        prime = next_prime
        @primes << prime
        @residues << modular_determinant(integer_matrix, prime)
        @modulus *= prime
        next unless @modulus > 2 * @bound

        combined = crt(@primes, @residues)
        return combined > @bound ? combined - @modulus : combined
      end
    end

    def next_prime
      prime =
        begin
          @prime_source.next
        rescue StopIteration
          raise InsufficientModuli, "the supplied primes are exhausted and the modulus " \
                                    "#{@modulus} does not exceed twice the bound #{@bound}"
        end
      unless prime.is_a?(Integer) && prime.prime?
        raise InvalidModulus, "#{prime} is not a prime number"
      end
      raise DuplicateModulus, "#{prime} was already used" if @primes.include?(prime)

      prime
    end

    # Gaussian elimination over the finite field GF(prime). A zero residue
    # is a valid result: it only means the prime divides the determinant.
    def modular_determinant(integer_matrix, prime)
      size = integer_matrix.size
      rows = integer_matrix.map { |row| row.map { |entry| entry % prime } }
      determinant = 1
      size.times do |column|
        pivot = (column...size).find { |row| rows[row][column] != 0 }
        return 0 if pivot.nil?

        if pivot != column
          rows[column], rows[pivot] = rows[pivot], rows[column]
          determinant = -determinant
        end
        determinant = determinant * rows[column][column] % prime
        eliminate_column(rows, column, prime)
      end
      determinant
    end

    def eliminate_column(rows, column, prime)
      inverse = mod_inverse(rows[column][column], prime)
      ((column + 1)...rows.size).each do |row|
        next if rows[row][column].zero?

        factor = rows[row][column] * inverse % prime
        (column...rows.size).each do |j|
          rows[row][j] = (rows[row][j] - factor * rows[column][j]) % prime
        end
      end
    end

    # Scales each row by the least common multiple of the denominators of
    # its entries, so that det(matrix) = det(integer_matrix) / scale with
    # both the clearing and the final rescaling exact.
    def clear_denominators(matrix)
      scale = 1
      integer_matrix = matrix.map do |row|
        row_lcm = row.map(&:denominator).inject(1) { |lcm, den| lcm.lcm(den) }
        scale *= row_lcm
        row.map { |entry| (entry * row_lcm).to_i }
      end
      [integer_matrix, scale]
    end

    # Explicit Hadamard bound: |det| <= product of the Euclidean row norms,
    # rounded up to an integer so the bound itself is an exact Integer.
    def hadamard_bound(integer_matrix)
      integer_matrix.inject(1) do |product, row|
        norm_squared = row.inject(0) { |sum, entry| sum + entry * entry }
        product * integer_sqrt_ceil(norm_squared)
      end
    end

    def integer_sqrt_ceil(number)
      root = Integer.sqrt(number)
      root * root == number ? root : root + 1
    end

    def validate_matrix(matrix)
      unless matrix.is_a?(Array) && !matrix.empty? && matrix.all? { |row| row.is_a?(Array) }
        raise InvalidMatrix, "the matrix must be a non-empty array of row arrays"
      end
      unless matrix.all? { |row| row.size == matrix.size }
        raise NonSquareMatrix, "expected a square matrix, got row sizes #{matrix.map(&:size)}"
      end
      matrix.each do |row|
        row.each do |entry|
          unless entry.is_a?(Integer) || entry.is_a?(Rational)
            raise InvalidMatrixEntry, "matrix entries must be Integer or Rational, got #{entry.class}"
          end
        end
      end
      matrix
    end
  end
end
