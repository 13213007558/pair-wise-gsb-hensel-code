# frozen_string_literal: true

require "test_helper"

class TestFormalPolynomial < Minitest::Test
  def test_initialize
    coefficients = [1, 2, 3]
    f = HenselCode::FormalPolynomial.new coefficients

    assert_equal HenselCode::FormalPolynomial, f.class
    assert_equal coefficients, f.coefficients
  end

  def test_initialize_trims_trailing_zeros
    f = HenselCode::FormalPolynomial.new [1, 2, 0, 0]

    assert_equal [1, 2], f.coefficients
  end

  def test_initialize_validations
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new nil }
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new 1 }
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new [] }
    assert_raises(ArgumentError) { HenselCode::FormalPolynomial.new [1, "x"] }
  end

  def test_zero_and_degree
    zero = HenselCode::FormalPolynomial.new [0]
    f = HenselCode::FormalPolynomial.new [1, 2, 3]

    assert_predicate zero, :zero?
    refute_predicate f, :zero?
    assert_equal 2, f.degree
  end

  def test_monic
    assert_predicate HenselCode::FormalPolynomial.new([2, 1]), :monic?
    refute_predicate HenselCode::FormalPolynomial.new([2, 3]), :monic?
    refute_predicate HenselCode::FormalPolynomial.new([0]), :monic?
  end

  def test_equality
    assert_equal HenselCode::FormalPolynomial.new([1, 2]), HenselCode::FormalPolynomial.new([1, 2, 0])
    refute_equal HenselCode::FormalPolynomial.new([1, 2]), HenselCode::FormalPolynomial.new([1, 3])
  end

  def test_addition_and_subtraction
    f = HenselCode::FormalPolynomial.new [1, 2, 3]
    g = HenselCode::FormalPolynomial.new [4, 5]

    assert_equal HenselCode::FormalPolynomial.new([5, 7, 3]), f + g
    assert_equal HenselCode::FormalPolynomial.new([-3, -3, 3]), f - g
  end

  def test_multiplication
    f = HenselCode::FormalPolynomial.new [1, 1]
    g = HenselCode::FormalPolynomial.new [2, 1]

    assert_equal HenselCode::FormalPolynomial.new([2, 3, 1]), f * g
    assert_equal HenselCode::FormalPolynomial.new([3, 3]), f * 3
  end

  def test_multiplication_with_invalid_type
    f = HenselCode::FormalPolynomial.new [1, 1]

    assert_raises(HenselCode::WrongHenselCodeInputType) { f * "x" }
  end

  def test_no_padic_carry
    # exact integer arithmetic: coefficients are not reduced or carried
    f = HenselCode::FormalPolynomial.new [9, 9]

    assert_equal HenselCode::FormalPolynomial.new([81, 162, 81]), f * f
  end

  def test_mod
    f = HenselCode::FormalPolynomial.new [7, -3, 12]

    assert_equal HenselCode::FormalPolynomial.new([2, 2, 2]), f.mod(5)
  end

  def test_divmod_mod
    f = HenselCode::FormalPolynomial.new [1, 0, 1] # x^2 + 1
    g = HenselCode::FormalPolynomial.new [1, 1]    # x + 1
    quotient, remainder = f.divmod_mod(g, 5)

    assert_equal HenselCode::FormalPolynomial.new([4, 1]), quotient
    assert_equal HenselCode::FormalPolynomial.new([2]), remainder
    assert_equal f.mod(5), ((quotient * g) + remainder).mod(5)
  end

  def test_gcd_mod
    f = HenselCode::FormalPolynomial.new [1, 2, 1] # (x+1)^2
    g = HenselCode::FormalPolynomial.new [1, 1]    # x+1

    assert_equal HenselCode::FormalPolynomial.new([1, 1]), f.gcd_mod(g, 5)
    assert_equal HenselCode::FormalPolynomial.new([1]), f.gcd_mod(HenselCode::FormalPolynomial.new([2, 1]), 5)
  end

  def test_inverse_mod
    g = HenselCode::FormalPolynomial.new [2, 1]      # x+2
    modulus = HenselCode::FormalPolynomial.new [1, 1, 1] # x^2+x+1, irreducible mod 5
    inverse = g.inverse_mod(modulus, 5)

    assert_equal HenselCode::FormalPolynomial.new([1]), (g * inverse).divmod_mod(modulus, 5).last
  end

  def test_inverse_mod_raises_when_not_coprime
    g = HenselCode::FormalPolynomial.new [1, 1]
    modulus = HenselCode::FormalPolynomial.new [1, 2, 1]

    assert_raises(ZeroDivisionError) { g.inverse_mod(modulus, 5) }
  end

  def test_to_s_and_inspect
    f = HenselCode::FormalPolynomial.new [2, 3, 4]

    assert_equal "2 + 3x + 4x^2", f.to_s
    assert_equal "<FormalPolynomial: 2 + 3x + 4x^2>", f.inspect
  end
end
