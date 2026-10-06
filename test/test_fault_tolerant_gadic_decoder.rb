# frozen_string_literal: true

require "test_helper"

class TestFaultTolerantGadicDecoder < Minitest::Test
  def setup
    @exponent = 2
    @primes = [257, 263, 269, 271, 277, 281]
    @moduli = @primes.map { |prime| prime**@exponent }
    @rational = Rational(37_897, 52_234)
    @numerator_bound = @rational.numerator
    @denominator_bound = @rational.denominator
    @residues = residues_for(@rational)
  end

  def residues_for(rational)
    @primes.map do |prime|
      HenselCode::TruncatedFinitePadicExpansion.new(prime, @exponent, rational).to_i
    end
  end

  def decode(residues, num_bound = @numerator_bound, den_bound = @denominator_bound, moduli = @moduli)
    HenselCode::FaultTolerantGadicDecoder.decode(moduli, residues, num_bound, den_bound)
  end

  def corrupt(residues, *indexes)
    corrupted = residues.dup
    indexes.each { |i| corrupted[i] = (corrupted[i] + 1) % @moduli[i] }
    corrupted
  end

  def test_recovers_rational_without_errors
    result = decode(@residues)

    assert_predicate result, :ok?
    assert_predicate result, :verified?
    assert_equal :ok, result.status
    assert_equal @rational, result.rational
    assert_equal @moduli, result.consistent_moduli
    assert_empty result.inconsistent_moduli
    assert_equal @moduli.inject(:*), result.consistent_product
    assert_equal 1, result.error_product
  end

  def test_recovers_rational_with_each_single_residue_corrupted
    @moduli.each_with_index do |modulus, i|
      result = decode(corrupt(@residues, i))

      assert_predicate result, :ok?, "failed to recover with residue #{i} corrupted"
      assert_equal @rational, result.rational
      assert_equal [modulus], result.inconsistent_moduli
      assert_equal @moduli.inject(:*) / modulus, result.consistent_product
      assert_equal modulus, result.error_product
    end
  end

  def test_recovers_rational_with_two_residues_corrupted
    rational = Rational(1234, 567)
    residues = residues_for(rational)

    (0...@moduli.size).to_a.combination(2) do |i, j|
      result = decode(corrupt(residues, i, j), rational.numerator, rational.denominator)

      assert_predicate result, :ok?, "failed to recover with residues #{i} and #{j} corrupted"
      assert_equal rational, result.rational
      assert_equal [@moduli[i], @moduli[j]], result.inconsistent_moduli
    end
  end

  def test_returns_checkable_consistent_modulus_product_condition
    result = decode(corrupt(@residues, 2))
    threshold = 2 * @numerator_bound * @denominator_bound * result.error_product

    assert_equal threshold, result.uniqueness_bound
    assert_operator result.consistent_product, :>, result.uniqueness_bound
    assert_equal result.consistent_product > result.uniqueness_bound, result.verified?
  end

  def test_does_not_accept_majority_consistent_rational_without_proof
    # 3/8 is consistent with 3 of the 4 residues (moduli 3, 5, 7), an outright
    # majority, and it is a small rational, yet the consistent-modulus product
    # 105 does not exceed 2 * 10 * 10 * 11, so it must not be certified.
    result = decode([0, 1, 3, 0], 10, 10, [3, 5, 7, 11])

    refute_predicate result, :ok?
    assert_equal :insufficient_information, result.status
    assert_nil result.rational
  end

  def test_reports_ambiguity_when_multiple_candidates_are_unprovable
    # residues are consistent with 1/2 on moduli {3, 5} and with 2/1 on moduli
    # {3, 7}; neither consistency set satisfies the product condition.
    result = decode([2, 3, 2], 2, 2, [3, 5, 7])

    assert_equal :ambiguous, result.status
    assert_nil result.rational
    assert_includes result.candidates, Rational(1, 2)
    assert_includes result.candidates, Rational(2, 1)
  end

  def test_reports_insufficient_information_when_unrecoverable
    # residues of 1/2 with the last two corrupted; the product of the two good
    # moduli (15) is too small to reconstruct or certify anything.
    result = decode([2, 3, 0, 0], 5, 5, [3, 5, 7, 11])

    refute_predicate result, :ok?
    assert_equal :insufficient_information, result.status
    assert_nil result.rational
  end

  def test_validates_inputs
    assert_raises(HenselCode::NonCoprimeModuli) { decode([1, 2], 1, 1, [4, 6]) }
    assert_raises(ArgumentError) { decode([1, 2], 1, 1, [3, 5, 7]) }
    assert_raises(ArgumentError) { decode([], 1, 1, []) }
    assert_raises(ArgumentError) { decode([3], 1, 1, [3, 5]) }
    assert_raises(ArgumentError) { decode([1, 2], 0, 1, [3, 5]) }
  end

  def test_truncated_finite_gadic_expansion_fault_tolerant_decode
    result = HenselCode::TruncatedFiniteGadicExpansion.fault_tolerant_decode(
      @moduli, corrupt(@residues, 4), @numerator_bound, @denominator_bound
    )

    assert_predicate result, :ok?
    assert_equal @rational, result.rational
    assert_equal [@moduli[4]], result.inconsistent_moduli
  end
end
