# frozen_string_literal: true

require "test_helper"

class TestMultiModularDeterminant < Minitest::Test
  # Exact fraction-free (Bareiss) reference determinant for integer matrices.
  def bareiss_determinant(matrix)
    size = matrix.size
    rows = matrix.map(&:dup)
    sign = 1
    previous_pivot = 1
    (0...size - 1).each do |column|
      pivot = (column...size).find { |row| rows[row][column] != 0 }
      return 0 if pivot.nil?

      if pivot != column
        rows[column], rows[pivot] = rows[pivot], rows[column]
        sign = -sign
      end
      ((column + 1)...size).each do |i|
        ((column + 1)...size).each do |j|
          rows[i][j] = (rows[i][j] * rows[column][column] - rows[i][column] * rows[column][j]) / previous_pivot
        end
      end
      previous_pivot = rows[column][column]
    end
    sign * rows[size - 1][size - 1]
  end

  def random_matrix(size, bits, seed)
    random = Random.new(seed)
    Array.new(size) do
      Array.new(size) { random.rand(-(2**bits)...2**bits) }
    end
  end

  def test_integer_matrix_matches_bareiss_reference
    matrix = random_matrix(5, 32, 42)
    assert_equal bareiss_determinant(matrix), HenselCode::MultiModularDeterminant.compute(matrix)
  end

  def test_large_integer_matrix_matches_bareiss_reference
    matrix = random_matrix(8, 64, 7)
    result = HenselCode::MultiModularDeterminant.new(matrix)
    assert_equal bareiss_determinant(matrix), result.determinant
    assert_operator result.modulus, :>, 2 * result.bound
    assert_equal result.primes.inject(:*), result.modulus
    assert_equal result.primes.size, result.residues.size
  end

  def test_negative_determinant
    assert_equal(-2, HenselCode::MultiModularDeterminant.compute([[1, 2], [3, 4]]))
  end

  def test_true_zero_determinant_is_not_confused_with_modular_singularity
    matrix = [
      [3, 1, 4],
      [1, 5, 9],
      [3, 1, 4]
    ]
    result = HenselCode::MultiModularDeterminant.new(matrix)
    assert_equal 0, result.determinant
    assert result.residues.all?(&:zero?), "every prime divides a zero determinant"
    assert_operator result.modulus, :>, 2 * result.bound
  end

  def test_zero_residues_under_first_primes_do_not_imply_singularity
    matrix = [
      [210, 0, 0],
      [0, 1, 0],
      [0, 0, 1]
    ]
    result = HenselCode::MultiModularDeterminant.new(matrix)
    assert_equal 210, result.determinant
    assert_equal [2, 3, 5, 7, 11], result.primes
    assert_equal [0, 0, 0, 0], result.residues.first(4)
    assert_equal 1, result.residues[4]
    assert_equal 210, result.bound
    assert_equal 2310, result.modulus
  end

  def test_bound_and_modulus_at_stopping_time
    result = HenselCode::MultiModularDeterminant.new([[3, 4], [1, 2]])
    assert_equal 2, result.determinant
    assert_equal 15, result.bound
    assert_equal [2, 3, 5, 7], result.primes
    assert_equal 210, result.modulus
  end

  def test_rational_input_with_exact_denominator_clearing
    matrix = [
      [Rational(1, 2), Rational(1, 3)],
      [Rational(1, 5), Rational(1, 7)]
    ]
    result = HenselCode::MultiModularDeterminant.new(matrix)
    assert_equal Rational(1, 210), result.determinant
    assert_equal 210, result.scale
    cleared = [[3, 2], [7, 5]]
    assert_equal bareiss_determinant(cleared), (result.determinant * result.scale).to_i
  end

  def test_rational_input_matching_integer_computation
    matrix = [
      [Rational(3, 2), Rational(-1, 4), 2],
      [Rational(1, 3), Rational(5, 6), Rational(-2, 3)],
      [1, Rational(-3, 5), Rational(7, 10)]
    ]
    result = HenselCode::MultiModularDeterminant.new(matrix)
    cleared = [
      [6, -1, 8],
      [2, 5, -4],
      [10, -6, 7]
    ]
    assert_equal 4 * 6 * 10, result.scale
    assert_equal Rational(bareiss_determinant(cleared), result.scale), result.determinant
  end

  def test_explicit_primes
    matrix = random_matrix(4, 8, 11)
    primes = [101, 103, 107, 109, 113, 127, 131]
    result = HenselCode::MultiModularDeterminant.new(matrix, primes: primes)
    assert_equal bareiss_determinant(matrix), result.determinant
    assert(result.primes.all? { |prime| primes.include?(prime) })
  end

  def test_class_level_compute
    assert_equal 1, HenselCode::MultiModularDeterminant.compute([[1]])
  end

  def test_invalid_matrix
    assert_raises(HenselCode::InvalidMatrix) { HenselCode::MultiModularDeterminant.compute([]) }
    assert_raises(HenselCode::InvalidMatrix) { HenselCode::MultiModularDeterminant.compute("nope") }
  end

  def test_non_square_matrix
    matrix = [[1, 2], [3, 4], [5, 6]]
    assert_raises(HenselCode::NonSquareMatrix) { HenselCode::MultiModularDeterminant.compute(matrix) }
  end

  def test_invalid_matrix_entry
    matrix = [[1, 2.5], [3, 4]]
    assert_raises(HenselCode::InvalidMatrixEntry) { HenselCode::MultiModularDeterminant.compute(matrix) }
  end

  def test_non_prime_modulus_is_rejected
    assert_raises(HenselCode::InvalidModulus) { HenselCode::MultiModularDeterminant.compute([[100]], primes: [4, 9]) }
  end

  def test_duplicate_prime_is_rejected
    assert_raises(HenselCode::DuplicateModulus) do
      HenselCode::MultiModularDeterminant.compute([[100]], primes: [2, 2, 3])
    end
  end

  def test_insufficient_primes_are_rejected
    error = assert_raises(HenselCode::InsufficientModuli) do
      HenselCode::MultiModularDeterminant.compute([[100]], primes: [2, 3])
    end
    assert_match(/exhausted/, error.message)
  end
end
