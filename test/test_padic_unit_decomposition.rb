# frozen_string_literal: true

require "test_helper"

class TestPadicUnitDecomposition < Minitest::Test
  def test_odd_prime_decomposition_reconstructs_input_and_lifts_root
    h = HenselCode::TruncatedFinitePadicExpansion.new(5, 4, 12)
    decomposition = h.unit_decomposition

    assert_equal :nonzero, decomposition[:status]
    assert_equal 0, decomposition[:valuation]
    assert_equal 625, decomposition[:unit_modulus]
    assert_equal 182, decomposition[:teichmuller_root]
    assert_equal 625, decomposition[:teichmuller_root_modulus]
    assert_equal 316, decomposition[:principal_unit]
    assert_equal 625, decomposition[:principal_unit_storage_modulus]
    refute_equal(12, decomposition[:teichmuller_root] % 5**4)
    assert_equal 1, decomposition[:principal_unit] % 5
    assert_equal 2, decomposition[:teichmuller_root] % 5
    assert_equal h.hensel_code, reconstruct(decomposition)
    assert_equal decomposition[:teichmuller_root], root_power(decomposition)
    assert_equal 1, root_order_power(decomposition)
    assert_equal :standard, decomposition[:teichmuller_convention]
    assert_equal :one_plus_padic_integers, decomposition[:principal_unit_convention]
  end

  def test_nonzero_valuation_preserves_effective_unit_precision
    h = HenselCode::TruncatedFinitePadicExpansion.new(5, 5, 300)
    decomposition = h.unit_decomposition

    assert_equal 2, decomposition[:valuation]
    assert_equal 25, decomposition[:valuation_power]
    assert_equal 3, decomposition[:unit_exponent]
    assert_equal 125, decomposition[:unit_modulus]
    assert_equal 12, decomposition[:unit]
    assert_equal 57, decomposition[:teichmuller_root]
    assert_equal 66, decomposition[:principal_unit]
    assert_equal 1, decomposition[:principal_unit] % 5
    assert_equal h.hensel_code, reconstruct(decomposition)
    assert_equal decomposition[:teichmuller_root], root_power(decomposition)
    assert_equal 1, root_order_power(decomposition)
  end

  def test_decompositions_are_compatible_when_input_precision_changes
    lower = HenselCode::TruncatedFinitePadicExpansion.new(5, 3, 12)
    higher = HenselCode::TruncatedFinitePadicExpansion.new(5, 5, 12)
    lower_decomposition = lower.unit_decomposition
    higher_decomposition = higher.unit_decomposition

    assert_equal lower_decomposition[:teichmuller_root] % lower_decomposition[:unit_modulus],
                 higher_decomposition[:teichmuller_root] % lower_decomposition[:unit_modulus]
    assert_equal lower_decomposition[:principal_unit] % lower_decomposition[:unit_modulus],
                 higher_decomposition[:principal_unit] % lower_decomposition[:unit_modulus]
    assert_equal 12, reconstruct(lower_decomposition)
    assert_equal 12, reconstruct(higher_decomposition)
  end

  def test_two_adic_signed_roots_and_four_principal_units
    positive = HenselCode::TruncatedFinitePadicExpansion.new(2, 4, 5)
    negative = HenselCode::TruncatedFinitePadicExpansion.new(2, 4, 14)
    low_precision = HenselCode::TruncatedFinitePadicExpansion.new(2, 1, 1)
    positive_decomposition = positive.unit_decomposition
    negative_decomposition = negative.unit_decomposition
    low_precision_decomposition = low_precision.unit_decomposition

    assert_equal 1, positive_decomposition[:teichmuller_root]
    assert_equal 5, positive_decomposition[:principal_unit]
    assert_equal 1, positive_decomposition[:principal_unit] % 4
    assert_equal :positive, positive_decomposition[:sign]
    assert_equal 1, root_power(positive_decomposition)
    assert_equal 1, root_order_power(positive_decomposition)
    assert_equal 5, reconstruct(positive_decomposition)

    assert_equal 1, negative_decomposition[:valuation]
    assert_equal 7, negative_decomposition[:teichmuller_root]
    assert_equal 1, negative_decomposition[:principal_unit]
    assert_equal :negative, negative_decomposition[:sign]
    assert_equal 1, root_power(negative_decomposition)
    assert_equal 1, root_order_power(negative_decomposition)
    assert_equal 14, reconstruct(negative_decomposition)

    assert_equal 1, low_precision_decomposition[:teichmuller_root]
    assert_equal :unsigned_low_precision, low_precision_decomposition[:teichmuller_convention]
    assert_equal :one_plus_two_adic_integers,
                 low_precision_decomposition[:principal_unit_convention]
    assert_equal 1, reconstruct(low_precision_decomposition)
  end

  def test_two_adic_decomposition_is_compatible_when_precision_changes
    lower = HenselCode::TruncatedFinitePadicExpansion.new(2, 2, 7)
    higher = HenselCode::TruncatedFinitePadicExpansion.new(2, 5, 7)
    lower_decomposition = lower.unit_decomposition
    higher_decomposition = higher.unit_decomposition

    assert_equal 3, lower_decomposition[:teichmuller_root]
    assert_equal 31, higher_decomposition[:teichmuller_root]
    assert_equal 1, lower_decomposition[:principal_unit]
    assert_equal 25, higher_decomposition[:principal_unit]
    assert_equal lower_decomposition[:teichmuller_root],
                 higher_decomposition[:teichmuller_root] % lower_decomposition[:unit_modulus]
    assert_equal lower_decomposition[:principal_unit],
                 higher_decomposition[:principal_unit] % lower_decomposition[:unit_modulus]
    assert_equal 3, reconstruct(lower_decomposition)
    assert_equal 7, reconstruct(higher_decomposition)
  end

  def test_exact_zero_and_residue_zero_have_separate_statuses
    exact_zero = HenselCode::TruncatedFinitePadicExpansion.new(5, 3, Rational(0, 1))
    residue_zero = HenselCode::TruncatedFinitePadicExpansion.new(5, 3, 125)

    assert_equal :exact_zero, exact_zero.unit_decomposition[:status]
    assert_equal :exact_zero, exact_zero.zero_status
    assert exact_zero.exact_zero?
    assert exact_zero.zero_modulo_prime_power?

    assert_equal :zero_modulo_prime_power, residue_zero.unit_decomposition[:status]
    assert_equal :zero_modulo_prime_power, residue_zero.zero_status
    refute residue_zero.exact_zero?
    assert residue_zero.zero_modulo_prime_power?
  end

  def test_integer_residue_inputs_are_normalized_without_reinterpreting_them
    h = HenselCode::TruncatedFinitePadicExpansion.new(5, 3, -1)

    assert_equal 124, h.hensel_code
    assert_equal :nonzero, h.unit_decomposition[:status]
    assert_equal 0, h.unit_decomposition[:valuation]
    assert_equal 124, reconstruct(h.unit_decomposition)
  end

  def test_finite_digit_decomposition_matches_residue_decomposition
    digit_expansion = HenselCode::FinitePadicExpansion.new(5, 3, [0, 2, 1])
    decomposition = digit_expansion.unit_decomposition

    assert_equal 35, decomposition[:residue]
    assert_equal 1, decomposition[:valuation]
    assert_equal 7, decomposition[:teichmuller_root]
    assert_equal 1, decomposition[:principal_unit]
    assert_equal [2, 1], decomposition[:teichmuller_root_digits]
    assert_equal [1, 0], decomposition[:principal_unit_digits]
    assert_equal 35, reconstruct(decomposition)
    assert_equal decomposition[:teichmuller_root], root_power(decomposition)
    assert_equal 1, root_order_power(decomposition)
  end

  private

  def reconstruct(decomposition)
    unit = (decomposition[:teichmuller_root] * decomposition[:principal_unit]) % decomposition[:unit_modulus]
    (decomposition[:valuation_power] * unit) % decomposition[:modulus]
  end

  def root_power(decomposition)
    modular_power(
      decomposition[:teichmuller_root],
      decomposition[:root_equation_exponent],
      decomposition[:unit_modulus]
    )
  end

  def root_order_power(decomposition)
    modular_power(
      decomposition[:teichmuller_root],
      decomposition[:root_order_exponent],
      decomposition[:unit_modulus]
    )
  end

  def modular_power(base, exponent, modulus)
    result = 1
    base %= modulus
    while exponent.positive?
      result = (result * base) % modulus if exponent.odd?
      base = (base * base) % modulus
      exponent /= 2
    end
    result
  end
end
