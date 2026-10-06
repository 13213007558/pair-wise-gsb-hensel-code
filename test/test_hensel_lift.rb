# frozen_string_literal: true

require "test_helper"

class TestHenselLift < Minitest::Test
  def setup
    @prime = 7
    @exponent = 5
    @modulus = @prime**@exponent
    # Coupled system: x and y cannot be solved for independently.
    # f1 = x^2 - y, f2 = y^2 - 2  =>  x^4 = 2 in the 7-adics, y = x^2.
    @equations = [
      HenselCode::MultivariatePolynomial.new([:x, :y], [2, 0] => 1, [0, 1] => -1),
      HenselCode::MultivariatePolynomial.new([:x, :y], [0, 2] => 1, [0, 0] => -2)
    ]
    @seed = { x: 2, y: 4 }
  end

  def lift(**overrides)
    HenselCode::HenselLift.lift(
      equations: overrides.fetch(:equations, @equations),
      prime: overrides.fetch(:prime, @prime),
      exponent: overrides.fetch(:exponent, @exponent),
      seed: overrides.fetch(:seed, @seed)
    )
  end

  def test_lift_coupled_system
    result = lift

    assert result.applicable?
    assert_equal :ok, result.status
    assert_equal @modulus, result.modulus
    assert result.verified?
  end

  def test_lifted_solution_satisfies_all_equations_modulo_p_to_the_n
    result = lift
    x = result.solution[:x]
    y = result.solution[:y]

    assert_equal 0, (x * x - y) % @modulus
    assert_equal 0, (y * y - 2) % @modulus
  end

  def test_lifted_solution_is_a_genuine_p_adic_solution
    result = lift
    x = result.solution[:x]
    y = result.solution[:y]

    # x^4 = 2 has no integer solution, so the lift must leave the seed behind
    refute_equal @seed[:x], x
    assert_equal 2, x**4 % @modulus
    # the coupling relation y = x^2 holds modulo p^N
    assert_equal y, x * x % @modulus
    # both unknowns reduce to the seed modulo p
    assert_equal @seed[:x], x % @prime
    assert_equal @seed[:y], y % @prime
  end

  def test_levels_report_residuals_and_effective_precision
    result = lift

    assert_equal @exponent, result.levels.size
    result.levels.each_with_index do |level, index|
      assert_equal index + 1, level[:precision]
      assert_equal @prime**(index + 1), level[:modulus]
      assert_equal [0, 0], level[:residuals]
    end
  end

  def test_final_residuals_are_substituted_back_into_all_equations
    result = lift

    assert_equal 2, result.residuals.size
    result.residuals.each do |residual|
      assert_equal 0, residual % @modulus
    end
  end

  def test_lift_with_integer_root
    # f1 = x*y - 6, f2 = x + y - 5 has the exact integer root (2, 3)
    equations = [
      HenselCode::MultivariatePolynomial.new([:x, :y], [1, 1] => 1, [0, 0] => -6),
      HenselCode::MultivariatePolynomial.new([:x, :y], [1, 0] => 1, [0, 1] => 1, [0, 0] => -5)
    ]
    result = lift(equations: equations, seed: { x: 2, y: 3 }, exponent: 3)

    assert result.applicable?
    assert_equal({ x: 2, y: 3 }, result.solution)
    assert result.verified?
  end

  def test_exponent_one_returns_the_seed
    result = lift(exponent: 1)

    assert result.applicable?
    assert_equal @seed, result.solution
    assert_equal 1, result.levels.size
    assert result.verified?
  end

  def test_singular_jacobian_is_reported_as_not_applicable
    # Jacobian of (x^2, y^2) at (0, 0) is the zero matrix modulo 7
    equations = [
      HenselCode::MultivariatePolynomial.new([:x, :y], [2, 0] => 1),
      HenselCode::MultivariatePolynomial.new([:x, :y], [0, 2] => 1)
    ]
    result = lift(equations: equations, seed: { x: 0, y: 0 })

    refute result.applicable?
    assert_equal :singular_jacobian, result.status
    assert_match(/not applicable/, result.message)
    assert_nil result.solution
    refute result.verified?
  end

  def test_seed_must_solve_the_system_modulo_p
    assert_raises(ArgumentError) { lift(seed: { x: 1, y: 1 }) }
  end

  def test_prime_must_be_prime
    assert_raises(ArgumentError) { lift(prime: 9) }
    assert_raises(ArgumentError) { lift(prime: 1) }
  end

  def test_exponent_must_be_positive
    assert_raises(ArgumentError) { lift(exponent: 0) }
  end

  def test_system_must_be_square
    assert_raises(ArgumentError) { lift(equations: [@equations.first]) }
  end

  def test_equations_must_share_variables
    other = HenselCode::MultivariatePolynomial.new([:x, :z], [1, 0] => 1)

    assert_raises(ArgumentError) { lift(equations: [@equations.first, other]) }
  end

  def test_seed_must_cover_all_variables
    assert_raises(ArgumentError) { lift(seed: { x: 2 }) }
  end
end
