# frozen_string_literal: true

require "test_helper"

class TestUnramifiedPadicExtension < Minitest::Test
  def setup
    @field = HenselCode::UnramifiedPadicExtension.new(2, 4, [1, 1, 1])
  end

  def test_degree_and_precision_are_independent
    field = HenselCode::UnramifiedPadicExtension.new(2, 7, [1, 1, 1])
    element = field.element([1, 1])

    assert_equal 2, field.degree
    assert_equal 7, field.precision
    assert_equal 2, element.coefficients.size
    assert_equal 7, element.precision
  end

  def test_addition_subtraction_multiplication_and_inverse
    left = @field.element([3, 5])
    right = @field.element([7, 11])

    assert_equal @field.element([10, 0]), left + right
    assert_equal @field.element([12, 10]), left - right
    assert_equal @field.element([1, 0]), left * left.inverse
    assert_equal left, left / right * right
  end

  def test_polynomial_reduction_and_coefficient_truncation_are_separate
    product = @field.element([1, 1]) * @field.element([1, 1])

    assert_equal [0, 1], product.coefficients
    assert_equal 4, product.precision

    lower_precision = HenselCode::UnramifiedPadicExtension.new(2, 3, [1, 1, 1]).element([1, 1])
    assert_equal [0, 7], lower_precision.frobenius.coefficients
  end

  def test_frobenius_returns_to_element_after_extension_degree
    units.each do |element|
      assert_equal element, element.frobenius(2)
      assert_equal element, element.frobenius(4)
      assert_equal element.frobenius, element.frobenius(3)
    end
  end

  def test_trace_and_norm_are_in_base_field_for_units
    units.each do |element|
      trace = element.trace
      norm = element.norm

      assert trace.base?
      assert norm.base?
      assert_equal element + element.frobenius, trace
      assert_equal element * element.frobenius, norm
    end
  end

  def test_quadratic_trace_and_norm_values
    element = @field.element([1, 1])

    assert_equal [1, 0], element.trace.coefficients
    assert_equal [1, 0], element.norm.coefficients
  end

  def test_elements_record_effective_precision
    unit = @field.element([1, 1])
    doubled = unit + unit
    restored = doubled - unit

    assert_equal 1, doubled.valuation
    assert_equal 3, doubled.coefficient_precision
    assert_equal 4, doubled.precision
    assert_equal unit, restored
  end

  def test_explicit_valuation_supports_p_inverse
    prime_element = @field.element([1, 0], valuation: 1, precision: 4)
    inverse = prime_element.inverse
    product = prime_element * inverse

    assert_equal 1, prime_element.valuation
    assert_equal(-1, inverse.valuation)
    assert_equal [1, 0], product.coefficients
    assert_equal 3, product.precision
  end

  def test_addition_respects_negative_valuation
    prime_inverse = @field.element([1, 0], valuation: -1, precision: 3)
    addition = prime_inverse + @field.one
    doubled = prime_inverse + prime_inverse

    assert_equal(-1, addition.valuation)
    assert_equal 4, addition.coefficient_precision
    assert_equal [3, 0], addition.coefficients
    assert_equal 0, doubled.valuation
    assert_equal 3, doubled.precision
  end

  def test_normalization_finds_multiple_powers_of_prime
    prime_squared = @field.element([4, 0])

    assert_equal 2, prime_squared.valuation
    assert_equal 2, prime_squared.coefficient_precision
    assert_equal 4, prime_squared.precision
  end

  def test_rejects_non_prime
    assert_raises(ArgumentError) { HenselCode::UnramifiedPadicExtension.new(4, 3, [1, 1, 1]) }
  end

  def test_rejects_non_positive_precision
    assert_raises(ArgumentError) { HenselCode::UnramifiedPadicExtension.new(2, 0, [1, 1, 1]) }
  end

  def test_rejects_non_monic_polynomial
    error = HenselCode::InvalidUnramifiedPolynomial
    assert_raises(error) { HenselCode::UnramifiedPadicExtension.new(2, 3, [1, 1, 2]) }
  end

  def test_rejects_reducible_polynomial
    error = HenselCode::InvalidUnramifiedPolynomial
    assert_raises(error) { HenselCode::UnramifiedPadicExtension.new(3, 3, [0, 1, 1]) }
    assert_raises(error) { HenselCode::UnramifiedPadicExtension.new(2, 3, [1, 0, 1]) }
  end

  def test_rejects_mixed_fields
    other_field = HenselCode::UnramifiedPadicExtension.new(3, 4, [1, 0, 1])

    assert_raises(HenselCode::IncompatibleExtensionFields) do
      @field.element([1, 1]) + other_field.element([1, 1])
    end
  end

  def test_rejects_zero_inverse
    assert_raises(ZeroDivisionError) { @field.zero.inverse }
  end

  private

  def units
    (0...2**4).each_with_object([]) do |constant, result|
      (0...2**4).each do |coefficient|
        element = @field.element([constant, coefficient])
        result << element unless element.zero? || !element.valuation.zero?
      end
    end
  end
end
