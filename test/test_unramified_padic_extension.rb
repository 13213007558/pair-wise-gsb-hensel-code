# frozen_string_literal: true

require "test_helper"

class TestUnramifiedPadicExtension < Minitest::Test
  def setup
    @prime = 5
    @precision = 6
    @modulus = [1, 1, 1] # theta^2 + theta + 1, irreducible modulo 5
    @field = HenselCode::UnramifiedPadicExtension.new(@prime, @modulus, @precision)
    @modulus_p = @prime**@precision
    @rng = Random.new(1234)
  end

  def random_element(precision: @precision)
    @field.element(Array.new(@field.degree) { @rng.rand(0...@modulus_p) }, precision: precision)
  end

  def test_initialize
    assert_equal @prime, @field.prime
    assert_equal @modulus, @field.modulus
    assert_equal @precision, @field.precision
    assert_equal 2, @field.degree
  end

  def test_rejects_non_prime
    [1, 4, 9, 15].each do |non_prime|
      assert_raises(ArgumentError) { HenselCode::UnramifiedPadicExtension.new(non_prime, @modulus, @precision) }
    end
  end

  def test_rejects_invalid_precision
    [0, -1, 1.5, "6"].each do |bad_precision|
      assert_raises(ArgumentError) { HenselCode::UnramifiedPadicExtension.new(@prime, @modulus, bad_precision) }
    end
  end

  def test_rejects_invalid_modulus_types
    assert_raises(ArgumentError) { HenselCode::UnramifiedPadicExtension.new(@prime, [1], @precision) }
    assert_raises(ArgumentError) { HenselCode::UnramifiedPadicExtension.new(@prime, [1, "1"], @precision) }
    assert_raises(ArgumentError) { HenselCode::UnramifiedPadicExtension.new(@prime, nil, @precision) }
  end

  def test_rejects_non_monic_modulus
    assert_raises(HenselCode::ModulusNotMonic) do
      HenselCode::UnramifiedPadicExtension.new(@prime, [1, 1, 2], @precision)
    end
  end

  def test_rejects_reducible_modulus_mod_p
    reducible_moduli = [
      [-1, 0, 1], # theta^2 - 1 = (theta - 1)(theta + 1)
      [0, 0, 1],  # theta^2
      [1, 2, 1]   # (theta + 1)^2
    ]
    reducible_moduli.each do |modulus|
      assert_raises(HenselCode::ModulusNotUnramified) do
        HenselCode::UnramifiedPadicExtension.new(@prime, modulus, @precision)
      end
    end
  end

  def test_accepts_degree_one_modulus
    field = HenselCode::UnramifiedPadicExtension.new(@prime, [3, 1], @precision)
    a = field.element([7])

    assert_equal 1, field.degree
    assert_equal a, a.frobenius
    assert_equal a, a.trace
    assert_equal a, a.norm
  end

  def test_basis_degree_is_not_padic_precision
    element = @field.element([1, 2, 3, 4, 5])

    assert_equal @field.degree, element.coefficients.size
    assert_equal @precision, element.precision
  end

  def test_modulus_reduction_and_truncation_are_separate
    # theta^2 = -1 - theta exactly; coefficients are truncated mod 5^6 afterwards
    theta_squared = @field.element([0, 0, 1])

    assert_equal [@modulus_p - 1, @modulus_p - 1], theta_squared.coefficients
    assert_equal [-1, -1], @field.reduce_modulus([0, 0, 1])
    assert_equal [@modulus_p - 1, @modulus_p - 1], @field.truncate([-1, -1], @precision)
  end

  def test_addition_and_subtraction
    a = random_element
    b = random_element

    assert_equal a, (a + b) - b
    assert_equal b, (a + b) - a
    assert @field.zero.zero?
    assert_equal a, a + @field.zero
  end

  def test_multiplication
    theta = @field.generator
    theta_squared = theta * theta

    assert_equal [@modulus_p - 1, @modulus_p - 1], theta_squared.coefficients
    assert_equal @field.one, theta * theta_squared # theta^3 = 1
  end

  def test_inverse
    a = random_element
    a = @field.element([1, 1]) until a.unit?

    assert_equal @field.one, a * a.inverse
    assert_equal a, a.inverse.inverse
  end

  def test_division
    a = random_element
    b = random_element
    b = @field.element([1, 1]) until b.unit?

    assert_equal a, (a / b) * b
  end

  def test_inverse_rejects_non_unit
    assert_raises(HenselCode::NonInvertibleElement) { @field.element([5, 10]).inverse }
    assert_raises(HenselCode::NonInvertibleElement) { @field.zero.inverse }
  end

  def test_frobenius_of_generator
    # sigma(theta) = theta^2 = -1 - theta exactly
    assert_equal [@modulus_p - 1, @modulus_p - 1], @field.generator.frobenius.coefficients
  end

  def test_frobenius_reduces_to_power_map_mod_p
    a = random_element
    sigma_a_mod_p = a.frobenius.coefficients.map { |c| c % @prime }
    # a^5 mod (modulus, 5) computed coefficient-wise over F_5
    power = @field.element(a.coefficients.map { |c| c % @prime }, precision: 1)
    4.times { power = power * @field.element(a.coefficients.map { |c| c % @prime }, precision: 1) }

    assert_equal power.coefficients.map { |c| c % @prime }, sigma_a_mod_p
  end

  def test_frobenius_repeated_degree_times_is_identity
    10.times do
      a = random_element
      b = a
      @field.degree.times { b = b.frobenius }

      assert_equal a, b
    end
  end

  def test_frobenius_is_ring_homomorphism
    a = random_element
    b = random_element

    assert_equal a.frobenius + b.frobenius, (a + b).frobenius
    assert_equal a.frobenius * b.frobenius, (a * b).frobenius
  end

  def test_trace_belongs_to_base_field
    10.times do
      trace = random_element.trace

      assert_equal 0, trace.coefficients[1] % @modulus_p
      assert_equal trace, trace.frobenius
    end
  end

  def test_norm_belongs_to_base_field
    10.times do
      norm = random_element.norm

      assert_equal 0, norm.coefficients[1] % @modulus_p
      assert_equal norm, norm.frobenius
    end
  end

  def test_trace_and_norm_exact_values
    # trace(a0 + a1*theta) = 2*a0 - a1 and norm(a0 + a1*theta) = a0^2 - a0*a1 + a1^2
    a = @field.element([3, 2])

    assert_equal [(2 * 3) - 2, 0], a.trace.coefficients
    assert_equal [(3**2) - (3 * 2) + (2**2), 0], a.norm.coefficients
    assert_equal [@modulus_p - 1, 0], @field.generator.trace.coefficients
    assert_equal [1, 0], @field.generator.norm.coefficients
  end

  def test_trace_additive_and_norm_multiplicative
    a = random_element
    b = random_element

    assert_equal a.trace + b.trace, (a + b).trace
    assert_equal a.norm * b.norm, (a * b).norm
  end

  def test_precision_tracking
    a = @field.element([1, 2], precision: 3)
    b = @field.element([3, 4])

    assert_equal 3, (a + b).precision
    assert_equal 3, (a - b).precision
    assert_equal 3, (a * b).precision
    assert_equal 3, (a / b).precision
    assert_equal 3, a.inverse.precision
    assert_equal 3, a.frobenius.precision
    assert_equal 3, a.trace.precision
    assert_equal 3, a.norm.precision
    assert_equal @precision, b.inverse.precision
  end

  def test_precision_tracking_coefficients
    a = @field.element([1, 2], precision: 3)
    b = @field.element([3, 4])
    product = a * b
    modulus_3 = @prime**3

    assert_equal ((1 * 3) - (2 * 4)) % modulus_3, product.coefficients[0]
    assert_equal ((1 * 4) + (2 * 3) - (2 * 4)) % modulus_3, product.coefficients[1]
  end

  def test_rejects_element_precision_above_field_precision
    assert_raises(ArgumentError) { @field.element([1, 2], precision: @precision + 1) }
  end

  def test_rejects_incompatible_operands
    other_field = HenselCode::UnramifiedPadicExtension.new(3, [1, 0, 1], @precision)
    a = @field.element([1, 2])
    b = other_field.element([1, 2])

    assert_raises(HenselCode::IncompatibleOperandTypes) { a + b }
    assert_raises(HenselCode::IncompatibleOperandTypes) { a * b }
  end

  def test_to_s_and_inspect
    a = @field.element([1, 2])

    assert_equal "1 + 2*theta (mod 5^6)", a.to_s
    assert_equal "<UnramifiedPadicExtension::Element: 1 + 2*theta (mod 5^6)>", a.inspect
  end
end
