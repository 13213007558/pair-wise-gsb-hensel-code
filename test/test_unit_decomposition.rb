# frozen_string_literal: true

require "test_helper"

class TestUnitDecomposition < Minitest::Test
  def test_truncated_odd_prime_decomposition_recombines_and_satisfies_equations
    h = HenselCode::TruncatedFinitePadicExpansion.new(3, 5, Rational(36, 7))
    decomposition = h.unit_decomposition

    assert_equal :nonzero, decomposition.status
    assert_equal 2, decomposition.valuation
    assert_equal 9, decomposition.valuation_power
    assert_equal 3, decomposition.unit_exponent
    assert_equal 1, decomposition.sign
    assert_equal :not_applicable, decomposition.sign_status
    assert_equal 1, decomposition.teichmuller_root.residue.pow(2, decomposition.unit_modulus)
    assert_equal 1, decomposition.principal_unit.residue % 3
    assert_equal h.to_i, decomposition.recombine.to_i
  end

  def test_finite_expansion_decomposition_uses_unit_precision
    h = HenselCode::FinitePadicExpansion.new(3, 5, Rational(36, 7))
    decomposition = h.unit_decomposition
    residue = 36 * HenselCode::TruncatedFinitePadicExpansion.new(3, 5, Rational(1, 7)).to_i
    residue %= 3**5

    assert_equal 2, decomposition.valuation
    assert_equal 3, decomposition.unit_exponent
    assert_equal 3, decomposition.teichmuller_root.exponent
    assert_equal 3, decomposition.principal_unit.exponent
    assert_equal residue, h.residue
    assert_equal h.to_a, decomposition.recombine.to_a
    assert_equal 1, decomposition.teichmuller_root.residue.pow(2, 27)
  end

  def test_teichmuller_root_is_lifted_beyond_first_digit
    h = HenselCode::TruncatedFinitePadicExpansion.new(3, 4, Rational(1, 5))
    root = h.unit_decomposition.teichmuller_root

    refute_equal 2, root.to_i
    assert_equal 80, root.to_i
    assert_equal 1, root.to_i.pow(2, 3**4)
  end

  def test_decompositions_are_compatible_when_input_precision_changes
    low = HenselCode::TruncatedFinitePadicExpansion.new(7, 3, Rational(11, 3))
    high = HenselCode::TruncatedFinitePadicExpansion.new(7, 5, Rational(11, 3))
    low_decomposition = low.unit_decomposition
    high_decomposition = high.unit_decomposition

    assert_equal low.to_i, low_decomposition.recombine.to_i
    assert_equal high.to_i, high_decomposition.recombine.to_i
    assert high_decomposition.compatible_with?(low_decomposition)
    assert low_decomposition.compatible_with?(high)
  end

  def test_p_two_uses_sign_and_four_principal_units
    h = HenselCode::TruncatedFinitePadicExpansion.new(2, 6, 56)
    decomposition = h.unit_decomposition

    assert_equal 3, decomposition.valuation
    assert_equal 3, decomposition.unit_exponent
    assert_equal(-1, decomposition.sign)
    assert_equal :determined, decomposition.sign_status
    assert_equal 7, decomposition.teichmuller_root.to_i
    assert_equal 1, decomposition.principal_unit.to_i
    assert_equal 1, decomposition.principal_unit.to_i % 4
    assert_equal 1, decomposition.teichmuller_root.to_i.pow(2, 8)
    assert_equal 56, decomposition.recombine.to_i
  end

  def test_p_two_decompositions_with_different_precision_are_compatible
    low = HenselCode::TruncatedFinitePadicExpansion.new(2, 4, Rational(9, 5))
    high = HenselCode::TruncatedFinitePadicExpansion.new(2, 6, Rational(9, 5))
    low_decomposition = low.unit_decomposition
    high_decomposition = high.unit_decomposition

    assert_equal 1, low_decomposition.sign
    assert_equal 1, high_decomposition.sign
    assert high_decomposition.compatible_with?(low_decomposition)
    assert_equal low.to_i, low_decomposition.recombine.to_i
    assert_equal high.to_i, high_decomposition.recombine.to_i
  end

  def test_exact_zero_and_zero_only_modulo_prime_power_have_distinct_statuses
    exact_zero = HenselCode::TruncatedFinitePadicExpansion.new(5, 3, Rational(0))
    residue_zero = HenselCode::TruncatedFinitePadicExpansion.new(5, 3, 0)
    prime_power_zero = HenselCode::TruncatedFinitePadicExpansion.new(5, 3, Rational(125))
    array_zero = HenselCode::FinitePadicExpansion.new(5, 3, [0, 0, 0])

    assert_equal :zero, exact_zero.zero_status
    assert exact_zero.zero?
    refute residue_zero.zero?
    assert_nil exact_zero.valuation
    assert_equal :zero, exact_zero.unit_decomposition.status
    assert_equal :zero_mod_prime_power, residue_zero.zero_status
    assert_nil residue_zero.valuation
    assert_equal :zero_mod_prime_power, residue_zero.unit_decomposition.status
    assert_equal :zero_mod_prime_power, prime_power_zero.zero_status
    assert_equal :zero_mod_prime_power, array_zero.zero_status
  end

  def test_p_two_one_digit_unit_does_not_invent_a_sign
    decomposition = HenselCode::TruncatedFinitePadicExpansion.new(2, 1, 1).unit_decomposition

    assert_nil decomposition.sign
    assert_equal :undetermined, decomposition.sign_status
    assert_equal 1, decomposition.teichmuller_root.to_i
    assert_equal 1, decomposition.principal_unit.to_i
    assert_equal 1, decomposition.recombine.to_i
  end

  def test_decomposition_hash_and_aliases_expose_components
    decomposition = HenselCode::TruncatedFinitePadicExpansion.new(3, 5, Rational(36, 7)).unit_decomposition

    assert_equal decomposition.teichmuller_root, decomposition[:teichmuller_root]
    assert_equal decomposition.teichmuller_root, decomposition.omega
    assert_equal decomposition.principal_unit, decomposition.principal
    assert_equal decomposition.valuation_power, decomposition.valuation_component
    assert_equal decomposition.status, decomposition.state
    assert_kind_of Hash, decomposition.to_h
  end
end
