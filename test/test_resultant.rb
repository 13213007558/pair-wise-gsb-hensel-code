# frozen_string_literal: true

require "test_helper"

class TestResultant < Minitest::Test
  def setup
    @prime = 257
  end

  def formal(coefficients)
    HenselCode::FormalPolynomial.new(@prime, coefficients)
  end

  def padic(integer, exponent)
    HenselCode::TruncatedFinitePadicExpansion.new(@prime, exponent, Rational(integer % (@prime**exponent), 1))
  end

  def padic_formal(coefficients, exponent)
    formal(coefficients.map { |coefficient| padic(coefficient, exponent) })
  end

  def test_resultant_of_linear_polynomials
    # res(x - a, x - b) = a - b
    assert_equal(-2, formal([-5, 1]).resultant(formal([-7, 1])).to_i)
  end

  def test_resultant_of_monic_polynomials
    # f = (x - 1)(x + 3), g monic cubic: res(f, g) = g(1) * g(-3)
    f = formal([-3, 2, 1])
    g = formal([5, 2, 0, 1])

    assert_equal 8 * (-28), f.resultant(g).to_i
  end

  def test_resultant_with_constant_polynomial
    # res(a, g) = a^deg(g)
    assert_equal 9, formal([3]).resultant(formal([1, 2, 1])).to_i
  end

  def test_resultant_antisymmetry
    f = formal([-3, 2, 1])
    g = formal([5, 2, 0, 1])
    h = formal([-2, 1])

    # deg(f) * deg(g) = 6 is even, deg(h) * deg(g) = 3 is odd
    assert_equal f.resultant(g).to_i, g.resultant(f).to_i
    assert_equal h.resultant(g).to_i, -g.resultant(h).to_i
  end

  def test_exact_resultant_is_determined
    result = formal([-3, 2, 1]).resultant(formal([5, 2, 0, 1]))

    assert result.exact?
    assert result.determined?
    assert_equal(-224, result.value)
    assert_equal(-224, result.to_i)
    assert_equal 0, result.valuation
    assert_equal 0, result.valuation_bound
    assert_equal Float::INFINITY, result.unit_precision
    refute result.zero?
    assert_equal false, result.common_factor?
  end

  def test_exact_common_factor
    # f = (x - 1)(x + 3), g = (x - 1)(x + 2)
    result = formal([-3, 2, 1]).resultant(formal([-2, 1, 1]))

    assert_equal 0, result.to_i
    assert result.zero?
    assert_equal Float::INFINITY, result.valuation
    assert_equal true, result.common_factor?
  end

  def test_padic_resultant_matches_integer_resultant_modulo_precision
    exponent = 5
    f = formal([-3, 2, 1])
    g = formal([5, 2, 0, 1])
    exact = f.resultant(g)
    result = padic_formal(f.coefficients, exponent).resultant(padic_formal(g.coefficients, exponent))

    refute result.exact?
    assert_equal exponent, result.precision
    assert_equal @prime**exponent, result.modulus
    assert_equal exact.to_i % (@prime**exponent), result.to_i
    assert result.determined?
    assert_equal false, result.common_factor?
  end

  def test_mixed_exact_and_padic_coefficients
    exponent = 4
    f = formal([padic(-3, exponent), 2, padic(1, exponent + 2)])
    g = formal([5, padic(2, exponent), 0, 1])
    exact = formal([-3, 2, 1]).resultant(formal([5, 2, 0, 1]))
    result = f.resultant(g)

    assert_equal exponent, result.precision
    assert_equal exact.to_i % (@prime**exponent), result.to_i
  end

  def test_vanishing_resultant_modulo_precision_is_not_a_common_factor
    # f and g share no common factor, but their resultant is divisible by
    # @prime**3, so modulo p^3 the resultant vanishes: no conclusion is allowed
    difference = @prime**3
    f = formal([-3, 2, 1])
    g = formal([-2 * (1 + difference), 1 - difference, 1])
    exact = f.resultant(g)

    assert_equal 3, exact.valuation

    result = padic_formal(f.coefficients, 3).resultant(padic_formal(g.coefficients, 3))

    assert_equal 0, result.to_i
    refute result.determined?
    refute result.zero?
    assert_nil result.valuation
    assert_equal 3, result.valuation_bound
    assert_equal 0, result.unit_precision
    assert_nil result.common_factor?
  end

  def test_genuine_common_factor_is_uncertain_at_finite_precision
    # even when a common factor really exists, finite precision can only
    # report that the resultant vanishes to the current precision
    f = formal([-3, 2, 1])
    g = formal([-2, 1, 1])
    result = padic_formal(f.coefficients, 4).resultant(padic_formal(g.coefficients, 4))

    assert_equal 0, result.to_i
    assert_nil result.common_factor?
    assert_nil result.valuation
    assert_equal 4, result.valuation_bound
  end

  def test_subresultant_of_degree_zero_equals_resultant
    10.times do
      f = Array.new(4) { rand(-20..20) }
      g = Array.new(3) { rand(-20..20) }
      next if f.last.zero? || g.last.zero?

      sequence = formal(f).subresultants(formal(g))

      assert_equal formal(f).resultant(formal(g)).to_i, sequence[0].coefficients.first
    end
  end

  def test_subresultants_detect_gcd_degree
    # gcd(f, g) = (x - 1)(x - 2) = x^2 - 3x + 2
    f = formal([-6, 11, -6, 1])
    g = formal([-8, 14, -7, 1])
    sequence = f.subresultants(g)

    assert_equal 4, sequence.size
    assert_equal [0, 1], sequence.vanishing_degrees
    assert_equal [2, 3], sequence.determined_degrees
    assert_equal 2, sequence.gcd_degree
    assert_equal [2, 2], sequence.gcd_degree_bounds

    leading = sequence[2].coefficients.last
    assert_equal 2 * leading, sequence[2].coefficients[0]
    assert_equal(-3 * leading, sequence[2].coefficients[1])
  end

  def test_subresultants_with_padic_coefficients
    exponent = 4
    f = [-6, 11, -6, 1]
    g = [-8, 14, -7, 1]
    exact_sequence = formal(f).subresultants(formal(g))
    sequence = padic_formal(f, exponent).subresultants(padic_formal(g, exponent))

    assert_equal exponent, sequence.precision
    assert_equal [0, 1], sequence.vanishing_degrees
    assert_equal [2, 3], sequence.determined_degrees
    assert_nil sequence.gcd_degree
    assert_equal [0, 2], sequence.gcd_degree_bounds

    s2 = sequence[2]
    exact_s2 = exact_sequence[2]
    s2.coefficients.each_with_index do |coefficient, index|
      assert_equal HenselCode::TruncatedFinitePadicExpansion, coefficient.class
      assert_equal exact_s2.coefficients[index] % (@prime**exponent), coefficient.hensel_code
    end
  end

  def test_subresultant_coefficient_valuations
    exponent = 4
    difference = @prime**2
    f = [-3, 2, 1]
    g = [-2 * (1 + difference), 1 - difference, 1]
    sequence = padic_formal(f, exponent).subresultants(padic_formal(g, exponent))

    assert_equal 2, sequence.coefficient_valuation(0, 0)
    assert_equal 2, sequence.valuation_bound(0)
  end

  def test_gcd_degree_zero_is_certified_by_nonzero_resultant
    difference = @prime**2
    f = [-3, 2, 1]
    g = [-2 * (1 + difference), 1 - difference, 1]
    sequence = padic_formal(f, 4).subresultants(padic_formal(g, 4))

    assert_equal 0, sequence.gcd_degree
  end

  def test_near_common_factor_stress_with_increasing_precision
    # f = (x - 1)(x + 3) and g = (x - 1 - p^8)(x + 2) share a factor only
    # modulo p^8; the exact resultant has valuation exactly 8
    difference_exponent = 8
    difference = @prime**difference_exponent
    f = [-3, 2, 1]
    g = [-2 * (1 + difference), 1 - difference, 1]
    exact = formal(f).resultant(formal(g))

    assert_equal difference_exponent, exact.valuation

    (1..12).each do |exponent|
      result = padic_formal(f, exponent).resultant(padic_formal(g, exponent))

      assert_equal exact.to_i % (@prime**exponent), result.to_i

      if exponent <= difference_exponent
        refute result.determined?
        assert_nil result.valuation
        assert_equal exponent, result.valuation_bound
        assert_nil result.common_factor?
      else
        assert result.determined?
        assert_equal difference_exponent, result.valuation
        assert_equal difference_exponent, result.valuation_bound
        assert_equal exponent - difference_exponent, result.unit_precision
        assert_equal false, result.common_factor?
      end
    end
  end
end
