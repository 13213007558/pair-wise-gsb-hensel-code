# frozen_string_literal: true

require "test_helper"

class TestMultiModular < Minitest::Test
  include HenselCode::MultiModular

  # Exact fraction-free Bareiss determinant, used as an independent reference.
  def reference_determinant(matrix)
    rows = matrix.map(&:dup)
    sign = 1
    previous = 1
    (0...rows.size - 1).each do |k|
      pivot = (k...rows.size).find { |i| rows[i][k] != 0 }
      return 0 if pivot.nil?

      unless pivot == k
        rows[k], rows[pivot] = rows[pivot], rows[k]
        sign = -sign
      end
      ((k + 1)...rows.size).each do |i|
        ((k + 1)...rows.size).each do |j|
          rows[i][j] = ((rows[i][j] * rows[k][k]) - (rows[i][k] * rows[k][j])) / previous
        end
      end
      previous = rows[k][k]
    end
    sign * rows[-1][-1]
  end

  def random_integer_matrix(size, bits, seed)
    rng = Random.new(seed)
    Array.new(size) { Array.new(size) { rng.rand(-(2**bits)...2**bits) } }
  end

  def test_hadamard_bound
    assert_equal 11, hadamard_bound([[1, 2], [3, 4]])
    assert_equal 0, hadamard_bound([[3, 4], [0, 0]])
    assert_equal 60, hadamard_bound([[3, 4], [12, 0]])
  end

  def test_determinant_mod
    matrix = random_integer_matrix(6, 32, 1234)
    expected = reference_determinant(matrix)
    [2, 3, 5, 257, 65_537].each do |prime|
      assert_equal expected % prime, determinant_mod(matrix, prime)
    end
  end

  def test_determinant_mod_singular_mod_prime_but_not_singular
    matrix = [[2, 0], [0, 1]]
    assert_equal 0, determinant_mod(matrix, 2)
    assert_equal 2, determinant_mod(matrix, 3)
  end

  def test_multi_modular_determinant_small
    matrix = [[1, 2], [3, 4]]
    result = multi_modular_determinant(matrix)
    assert_equal(-2, result.determinant)
  end

  def test_multi_modular_determinant_big_integer_entries
    matrix = random_integer_matrix(8, 64, 42)
    expected = reference_determinant(matrix)
    result = multi_modular_determinant(matrix)
    assert_equal expected, result.determinant
  end

  def test_multi_modular_determinant_negative
    matrix = random_integer_matrix(7, 48, 7)
    expected = reference_determinant(matrix)
    result = multi_modular_determinant(matrix)
    assert_equal expected, result.determinant
    assert result.modulus > 2 * result.bound
    assert result.bound >= expected.abs
  end

  def test_multi_modular_determinant_true_zero
    matrix = random_integer_matrix(6, 32, 99)
    matrix[3] = matrix[1].dup
    result = multi_modular_determinant(matrix)
    assert_equal 0, reference_determinant(matrix)
    assert_equal 0, result.determinant
    assert result.remainders.all?(&:zero?)
  end

  def test_multi_modular_determinant_zero_mod_first_primes_only
    matrix = [[2, 0, 0, 0, 0, 0],
              [0, 3, 0, 0, 0, 0],
              [0, 0, 5, 0, 0, 0],
              [0, 0, 0, 7, 0, 0],
              [0, 0, 0, 0, 11, 0],
              [0, 0, 0, 0, 0, 13]]
    result = multi_modular_determinant(matrix, primes: Prime.first(8))
    assert_equal 30_030, result.determinant
    assert_equal [0, 0, 0, 0, 0, 0], result.remainders.first(6)
    assert_equal 30_030, result.bound
    assert_equal 510_510, result.modulus
    assert result.modulus > 2 * result.bound
  end

  def test_multi_modular_determinant_rational
    matrix = [[Rational(1, 2), Rational(1, 3)], [Rational(1, 5), Rational(1, 7)]]
    result = multi_modular_determinant(matrix)
    assert_equal Rational(1, 210), result.determinant
  end

  def test_multi_modular_determinant_rational_random
    rng = Random.new(2024)
    matrix = Array.new(4) { Array.new(4) { Rational(rng.rand(-50..50), rng.rand(1..20)) } }
    expected = reference_determinant(matrix)
    result = multi_modular_determinant(matrix)
    assert_equal expected, result.determinant
  end

  def test_multi_modular_determinant_mixed_integer_and_rational
    matrix = [[2, Rational(1, 2)], [Rational(3, 4), 5]]
    result = multi_modular_determinant(matrix)
    assert_equal Rational(77, 8), result.determinant
  end

  def test_non_square_matrix_raises
    assert_raises(HenselCode::NonSquareMatrix) { multi_modular_determinant([[1, 2, 3], [4, 5, 6]]) }
    assert_raises(HenselCode::NonSquareMatrix) { multi_modular_determinant([]) }
    assert_raises(HenselCode::NonSquareMatrix) { determinant_mod([[1, 2, 3], [4, 5, 6]], 7) }
    assert_raises(HenselCode::NonSquareMatrix) { hadamard_bound([[1, 2, 3], [4, 5, 6]]) }
  end
end
