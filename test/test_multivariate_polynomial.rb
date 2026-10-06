# frozen_string_literal: true

require "test_helper"

class TestMultivariatePolynomial < Minitest::Test
  def setup
    @variables = [:x, :y]
    @f = HenselCode::MultivariatePolynomial.new(@variables, [2, 0] => 1, [0, 1] => -1) # x^2 - y
  end

  def test_initialize
    assert_equal [:x, :y], @f.variables
    assert_equal({ [2, 0] => 1, [0, 1] => -1 }, @f.terms)
  end

  def test_initialize_drops_zero_coefficients
    f = HenselCode::MultivariatePolynomial.new(@variables, [1, 0] => 2, [0, 1] => 0)

    assert_equal({ [1, 0] => 2 }, f.terms)
  end

  def test_initialize_validations
    assert_raises(ArgumentError) { HenselCode::MultivariatePolynomial.new([], [1, 0] => 1) }
    assert_raises(ArgumentError) { HenselCode::MultivariatePolynomial.new([:x, :x], [1, 0] => 1) }
    assert_raises(ArgumentError) { HenselCode::MultivariatePolynomial.new(@variables, nil) }
    assert_raises(ArgumentError) { HenselCode::MultivariatePolynomial.new(@variables, [1] => 1) }
    assert_raises(ArgumentError) { HenselCode::MultivariatePolynomial.new(@variables, [-1, 0] => 1) }
    assert_raises(ArgumentError) { HenselCode::MultivariatePolynomial.new(@variables, [1, 0] => 1.5) }
  end

  def test_evaluate
    assert_equal 0, @f.evaluate(x: 2, y: 4)
    assert_equal 5, @f.evaluate(x: 3, y: 4)
  end

  def test_evaluate_requires_full_point
    assert_raises(ArgumentError) { @f.evaluate(x: 2) }
    assert_raises(ArgumentError) { @f.evaluate(nil) }
  end

  def test_partial
    dfdx = @f.partial(:x)
    dfdy = @f.partial(:y)

    assert_equal HenselCode::MultivariatePolynomial.new(@variables, [1, 0] => 2), dfdx
    assert_equal HenselCode::MultivariatePolynomial.new(@variables, [0, 0] => -1), dfdy
  end

  def test_partial_of_constant_term_vanishes
    f = HenselCode::MultivariatePolynomial.new(@variables, [0, 0] => 7)

    assert f.partial(:x).zero?
  end

  def test_partial_unknown_variable
    assert_raises(ArgumentError) { @f.partial(:z) }
  end

  def test_degree
    assert_equal 2, @f.degree
    assert_equal 0, HenselCode::MultivariatePolynomial.new(@variables, {}).degree
  end

  def test_zero
    assert HenselCode::MultivariatePolynomial.new(@variables, {}).zero?
    refute @f.zero?
  end

  def test_equality
    assert_equal HenselCode::MultivariatePolynomial.new(@variables, [2, 0] => 1, [0, 1] => -1), @f
    refute_equal HenselCode::MultivariatePolynomial.new(@variables, [2, 0] => 1), @f
  end

  def test_to_s
    assert_equal "x^2 - y", @f.to_s

    g = HenselCode::MultivariatePolynomial.new(@variables, [1, 0] => 2, [0, 1] => 3, [0, 0] => -5)
    assert_equal "2*x + 3*y - 5", g.to_s
    assert_equal "0", HenselCode::MultivariatePolynomial.new(@variables, {}).to_s
  end

  def test_inspect
    assert_equal "<MultivariatePolynomial: x^2 - y>", @f.inspect
  end
end
