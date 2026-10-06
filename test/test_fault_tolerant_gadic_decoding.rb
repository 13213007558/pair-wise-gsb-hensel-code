# frozen_string_literal: true

require "test_helper"

class TestFaultTolerantGadicDecoding < Minitest::Test
  def setup
    @moduli = [7, 11, 13, 17]
    @numerator_bound = 10
    @denominator_bound = 10
    @rational = Rational(2, 3)
    @remainders = encode_remainders(@rational, @moduli)
  end

  def test_recovers_rational_when_one_remainder_is_corrupted
    @moduli.each_index do |corrupted_index|
      remainders = @remainders.dup
      remainders[corrupted_index] = (@remainders[corrupted_index] + 1) % @moduli[corrupted_index]

      result = decode(remainders)

      assert_equal :recovered, result[:status]
      assert_equal @rational, result[:rational]
      assert_equal [corrupted_index], result[:inconsistent_indices]
      assert_equal [@moduli[corrupted_index]], result[:inconsistent_moduli]
      expected_moduli = @moduli.reject.with_index { |_, index| index == corrupted_index }
      assert_equal expected_moduli, result[:consistent_moduli]
      assert_equal 200, result[:required_modulus_product]
      assert_operator result[:consistent_modulus_product], :>, result[:required_modulus_product]
    end
  end

  def test_recovers_rational_when_two_remainders_are_corrupted
    corrupted_indices = [0, 3]
    remainders = @remainders.dup
    remainders[0] = 0
    remainders[3] = 1

    result = decode(remainders, 2, 3)

    assert_equal :recovered, result[:status]
    assert_equal @rational, result[:rational]
    assert_equal corrupted_indices, result[:inconsistent_indices]
    assert_equal [11, 13], result[:consistent_moduli]
    assert_equal 143, result[:consistent_modulus_product]
    assert_equal 12, result[:required_modulus_product]
    assert_equal "143 > 12", result[:consistent_modulus_product_condition]
  end

  def test_reports_ambiguity_when_two_candidates_have_enough_consistent_moduli
    moduli = [7, 11, 13, 17, 19]
    remainders = [3, 8, 5, 7, 4]

    result = HenselCode::Tools.fault_tolerant_gadic_decode(moduli, remainders, 2, 5, max_errors: 2)
    rationals = result[:candidates].map { |candidate| candidate[:rational] }

    assert_equal :ambiguous, result[:status]
    assert_includes rationals, Rational(2, 3)
    assert_includes rationals, Rational(1, 5)
    assert_equal 20, result[:required_modulus_product]
    assert(result[:candidates].all? { |candidate| candidate[:consistent_modulus_product] > 20 })
  end

  def test_reports_unrecoverable_when_no_consistent_set_is_large_enough
    moduli = [3, 5, 7, 11]
    remainders = [0, 0, 4, 6]

    result = HenselCode::Tools.fault_tolerant_gadic_decode(moduli, remainders, 20, 20)

    assert_equal :unrecoverable, result[:status]
    refute result.key?(:rational)
    assert_equal 800, result[:required_modulus_product]
  end

  def test_accepts_separate_numerator_and_denominator_bounds
    result = decode(@remainders, 5, 3)

    assert_equal :recovered, result[:status]
    assert_equal @rational, result[:rational]
    assert_equal 5, result[:numerator_bound]
    assert_equal 3, result[:denominator_bound]
    assert_equal 30, result[:required_modulus_product]
  end

  def test_recovers_zero_with_an_alias
    result = HenselCode::Tools.fault_tolerant_decode(@moduli, [0, 0, 0, 1], 1, 1)

    assert_equal :recovered, result[:status]
    assert_equal Rational(0, 1), result[:rational]
    assert_equal [17], result[:inconsistent_moduli]
  end

  def test_validates_pairwise_coprime_moduli_and_remainder_range
    assert_raises(ArgumentError) { decode(@remainders, 10, 10, [7, 11, 13, 26]) }
    bad_remainders = @remainders.dup
    bad_remainders[0] = @moduli[0]

    assert_raises(ArgumentError) { decode(bad_remainders) }
  end

  private

  def decode(
    remainders,
    numerator_bound = @numerator_bound,
    denominator_bound = @denominator_bound,
    moduli = @moduli
  )
    HenselCode::Tools.fault_tolerant_gadic_decode(
      moduli,
      remainders,
      numerator_bound,
      denominator_bound
    )
  end

  def encode_remainders(rational, moduli)
    moduli.map do |modulus|
      (rational.numerator * HenselCode::Tools.mod_inverse(rational.denominator, modulus)) % modulus
    end
  end
end
