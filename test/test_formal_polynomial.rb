# frozen_string_literal: true

require "test_helper"

class TestFormalPolynomial < Minitest::Test
  def test_carry_free_arithmetic
    f = HenselCode::FormalPolynomial.new([3, 3])
    g = HenselCode::FormalPolynomial.new([3, 3])

    # no p-adic carry: coefficients may grow beyond any modulus
    assert_equal [9, 18, 9], (f * g).coefficients
    assert_equal [6, 6], (f + g).coefficients
    assert_equal [0], (f - g).coefficients
  end

  def test_coefficients_are_trimmed
    assert_equal [1, 2], HenselCode::FormalPolynomial.new([1, 2, 0, 0]).coefficients
  end

  def test_degree_and_monic
    f = HenselCode::FormalPolynomial.new([-6, 11, -6, 1])

    assert_equal 3, f.degree
    assert f.monic?
    refute HenselCode::FormalPolynomial.new([1, 2]).monic?
  end

  def test_divmod_by_monic_divisor
    f = HenselCode::FormalPolynomial.new([-6, 11, -6, 1])
    g = HenselCode::FormalPolynomial.new([-1, 1])
    quotient, remainder = f.divmod(g)

    assert_equal [6, -5, 1], quotient.coefficients
    assert remainder.zero?
  end

  def test_divmod_with_remainder
    f = HenselCode::FormalPolynomial.new([1, 0, 1])
    g = HenselCode::FormalPolynomial.new([-1, 1])
    quotient, remainder = f.divmod(g)

    assert_equal [1, 1], quotient.coefficients
    assert_equal [2], remainder.coefficients
  end

  def test_divmod_inexact_leading_coefficient_raises
    f = HenselCode::FormalPolynomial.new([1, 1])
    g = HenselCode::FormalPolynomial.new([1, 2])

    assert_raises(ArgumentError) { f.divmod(g) }
  end

  def test_divmod_mod
    f = HenselCode::FormalPolynomial.new([-6, 11, -6, 1])
    g = HenselCode::FormalPolynomial.new([-1, 1])
    quotient, remainder = f.divmod_mod(5, g)

    assert_equal [1, 0, 1], quotient.coefficients
    assert remainder.zero?
  end

  def test_gcd_mod
    f = HenselCode::FormalPolynomial.new([-1, 0, 1])
    g = HenselCode::FormalPolynomial.new([-1, 1])

    assert_equal HenselCode::FormalPolynomial.new([4, 1]), f.gcd_mod(5, g)
    assert_equal HenselCode::FormalPolynomial.one, f.gcd_mod(5, HenselCode::FormalPolynomial.new([-2, 1]))
  end

  def test_extended_gcd_mod
    a = HenselCode::FormalPolynomial.new([1, 0, 1])
    b = HenselCode::FormalPolynomial.new([1, 1])
    gcd, s, t = a.extended_gcd_mod(5, b)

    assert gcd.monic?
    assert (((s * a) + (t * b)) - gcd).mod(5).zero?
  end

  def test_inverse_mod
    a = HenselCode::FormalPolynomial.new([1, 1])
    modulus_polynomial = HenselCode::FormalPolynomial.new([1, 0, 1])
    inverse = a.inverse_mod(modulus_polynomial, 5)

    assert_equal HenselCode::FormalPolynomial.new([3, 2]), inverse
    assert_equal HenselCode::FormalPolynomial.one, (a * inverse).divmod_mod(5, modulus_polynomial).last
  end

  def test_inverse_mod_non_coprime_raises
    a = HenselCode::FormalPolynomial.new([2, 1])
    modulus_polynomial = HenselCode::FormalPolynomial.new([1, 0, 1])

    assert_raises(ArgumentError) { a.inverse_mod(modulus_polynomial, 5) }
  end

  def test_invalid_coefficients_raise
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new(nil) }
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new("1, 2") }
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new([1, "2"]) }
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new([]) }
  end
end
