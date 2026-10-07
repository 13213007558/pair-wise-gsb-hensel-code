# frozen_string_literal: true

require "test_helper"

class TestPolynomialResultant < Minitest::Test
  def test_formal_polynomial_is_independent_from_digit_carry_polynomial
    f = HenselCode::FormalPolynomial.new(5, 2, [1, 2])
    g = HenselCode::FormalPolynomial.new(5, 2, [3, 4])

    assert_equal [3, 10, 8].map { |coefficient| coefficient % 25 }, (f * g).coefficients
  end

  def test_resultant_over_integers
    f = HenselCode::FormalPolynomial.new(257, 3, [4, 3, 2])
    g = HenselCode::FormalPolynomial.new(257, 3, [6, 5])

    result = HenselCode::PolynomialResultant.new(f, g).resultant

    assert_equal 82, result[:residue]
    assert_equal 0, result[:valuation]
    assert_equal 0, result[:valuation_lower_bound]
    assert result[:valuation_determined]
    assert_equal 3, result[:known_unit_precision]
  end

  def test_subresultants_match_integer_subresultant_chain
    f = HenselCode::FormalPolynomial.new(257, 3, [1, 1, 0, 3])
    g = HenselCode::FormalPolynomial.new(257, 3, [4, -1, 2])

    subresultants = HenselCode::PolynomialResultant.new(f, g).subresultants
    first = subresultants[0][:coefficients].map { |coefficient| coefficient[:residue] }
    second = subresultants[1][:coefficients].map { |coefficient| coefficient[:residue] }

    assert_equal [355], first
    assert_equal [-8, -17].map { |coefficient| coefficient % 257**3 }, second
  end

  def test_mod_p_zero_does_not_claim_a_common_factor
    f = HenselCode::FormalPolynomial.from_padic_coefficients(5, 1, [[1], [1]])
    g = HenselCode::FormalPolynomial.from_padic_coefficients(5, 1, [[1], [1]])

    result = HenselCode::PolynomialResultant.new(f, g).resultant

    assert_equal 0, result[:residue]
    assert_nil result[:valuation]
    assert_equal 1, result[:valuation_lower_bound]
    refute result[:valuation_determined]
    assert_equal 0, result[:known_unit_precision]
  end

  def test_near_common_factor_recovers_integer_resultant_as_precision_grows
    exact_resultant = -(5**3)
    observed = (1..5).map do |precision|
      f = HenselCode::FormalPolynomial.new(5, precision, [1, 1])
      g = HenselCode::FormalPolynomial.new(5, precision, [1 + 5**3, 1])
      HenselCode::PolynomialResultant.new(f, g).resultant
    end

    observed.first(3).each do |result|
      assert_equal 0, result[:residue]
      assert_nil result[:valuation]
      refute result[:valuation_determined]
      assert_equal 0, result[:known_unit_precision]
      assert_equal result[:precision], result[:unresolved_precision]
    end

    determined = observed[3]
    assert_equal exact_resultant % (5**4), determined[:residue]
    assert_equal 3, determined[:valuation]
    assert_equal 3, determined[:valuation_lower_bound]
    assert determined[:valuation_determined]
    assert_equal 1, determined[:known_unit_precision]
    assert_equal 0, determined[:unresolved_precision]

    refined = observed.last
    assert_equal exact_resultant % (5**5), refined[:residue]
    assert_equal 3, refined[:valuation]
    assert_equal 2, refined[:known_unit_precision]
    assert_equal 0, refined[:unresolved_precision]
  end

  def test_near_common_factor_accepts_hensel_code_coefficients
    f_code = HenselCode::FinitePadicExpansion.new(5, 4, Rational(1, 1))
    x_code = HenselCode::FinitePadicExpansion.new(5, 4, Rational(1, 1))
    f = HenselCode::FormalPolynomial.from_padic_coefficients(5, 4, [f_code, x_code])
    g = HenselCode::FormalPolynomial.new(5, 4, [1 + 5**3, 1])

    result = HenselCode::PolynomialResultant.resultant(f, g)

    assert_equal(-(5**3) % 5**4, result[:residue])
    assert_equal 3, result[:valuation]
  end

  def test_subresultant_reports_coefficient_level_precision
    f = HenselCode::FormalPolynomial.new(5, 3, [1, 1])
    g = HenselCode::FormalPolynomial.new(5, 3, [1 + 5**3, 1])

    constant = HenselCode::PolynomialResultant.new(f, g).subresultants.first[:coefficients].first

    assert_equal 0, constant[:residue]
    assert_nil constant[:valuation]
    assert_equal 3, constant[:valuation_lower_bound]
    refute constant[:valuation_determined]
    assert_equal 0, constant[:known_unit_precision]
  end
end
