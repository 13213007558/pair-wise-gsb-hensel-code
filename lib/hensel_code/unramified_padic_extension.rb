# frozen_string_literal: true

module HenselCode
  # Unramified extension Q_p(alpha), where alpha is a root of an
  # irreducible monic polynomial over F_p.
  class UnramifiedPadicExtension
    include ModularPolynomialArithmetic

    attr_reader :prime, :precision, :degree, :polynomial

    def initialize(prime, precision, irreducible_polynomial)
      @polynomial = polynomial_coefficients(irreducible_polynomial)
      validate_arguments(prime, precision)
      @prime = prime
      @precision = precision
      validate_polynomial
      @degree = @polynomial.size - 1
      @frobenius_roots = {}
    end

    def element(coefficients, valuation: 0, precision: @precision)
      UnramifiedPadicExtensionElement.new(self, coefficients, valuation, precision)
    end

    alias [] element

    def zero
      build_element(Array.new(degree, 0), 0, precision, true)
    end

    def one
      element([1] + Array.new(degree - 1, 0))
    end

    def build_element(coefficients, valuation, coefficient_precision, zero = false)
      UnramifiedPadicExtensionElement.allocate.tap do |field_element|
        field_element.send(:initialize_internal, self, coefficients, valuation, coefficient_precision, zero)
      end
    end

    def reduction_polynomial(precision = @precision)
      reduce_coefficients(polynomial, prime**precision)
    end

    def reduce(coefficients, precision = @precision)
      modulus = prime**precision
      reduced = modulo(reduce_coefficients(coefficients, modulus), reduction_polynomial(precision), modulus)
      pad_to_degree(reduced.map { |coefficient| coefficient % modulus })
    end

    def multiply_unit(left, right, precision)
      multiply_mod(left, right, reduction_polynomial(precision), prime**precision)
    end

    def inverse_unit(coefficients, precision)
      modulus = prime**precision
      reduced = reduce(coefficients, precision)
      inverse = inverse_mod(reduced.map { |coefficient| coefficient % prime }, polynomial, prime)
      current_precision = 1

      while current_precision < precision
        next_precision = [current_precision * 2, precision].min
        next_modulus = prime**next_precision
        next_reduced = reduced.map { |coefficient| coefficient % next_modulus }
        next_inverse = inverse.map { |coefficient| coefficient % next_modulus }
        product = multiply_unit(next_inverse, next_reduced, next_precision)
        correction = multiply_unit(next_inverse, subtract([1], product, next_modulus), next_precision)
        inverse = add(next_inverse, correction, next_modulus)
        current_precision = next_precision
      end

      pad_to_degree(inverse.map { |coefficient| coefficient % modulus })
    end

    def frobenius_root(precision)
      @frobenius_roots[precision] ||= lift_frobenius_root(precision)
    end

    def same_configuration?(other)
      other.is_a?(UnramifiedPadicExtension) &&
        prime == other.prime &&
        precision == other.precision &&
        polynomial == other.polynomial
    end

    def evaluate(polynomial_coefficients, point, precision)
      modulus = prime**precision
      result = [0]
      polynomial_coefficients.reverse_each do |coefficient|
        result = multiply_unit(result, point, precision)
        result = add(result, [coefficient % modulus], modulus)
      end
      pad_to_degree(result.map { |coefficient| coefficient % modulus })
    end

    def pad_to_degree(coefficients)
      coefficients.dup.fill(0, coefficients.size...degree)
    end

    private

    def polynomial_coefficients(value)
      if value.is_a?(Polynomial)
        value.coefficients.dup
      elsif value.is_a?(Array)
        value.dup
      else
        raise InvalidUnramifiedPolynomial, "polynomial must be an Array or HenselCode::Polynomial"
      end
    end

    def validate_arguments(prime, precision)
      unless prime.is_a?(Integer) && Prime.prime?(prime)
        raise ArgumentError, "prime must be a prime integer"
      end
      unless precision.is_a?(Integer) && precision.positive?
        raise ArgumentError, "precision must be a positive integer"
      end
    end

    def validate_polynomial
      valid = @polynomial.is_a?(Array) && !@polynomial.empty? && @polynomial.all? { |c| c.is_a?(Integer) }
      raise InvalidUnramifiedPolynomial, "polynomial coefficients must be integers" unless valid

      @polynomial = trim(@polynomial)
      if @polynomial.size <= 1
        raise InvalidUnramifiedPolynomial, "polynomial must have positive degree"
      end
      unless @polynomial.last % prime == 1
        raise InvalidUnramifiedPolynomial, "polynomial must be monic"
      end
      unless irreducible_over_prime_field?(@polynomial, prime)
        raise InvalidUnramifiedPolynomial, "polynomial must be irreducible modulo the prime"
      end
    end

    def lift_frobenius_root(target_precision)
      root = power_mod([0, 1], prime, polynomial, prime)
      current_precision = 1

      while current_precision < target_precision
        next_precision = [current_precision * 2, target_precision].min
        modulus = prime**next_precision
        polynomial = reduction_polynomial(next_precision)
        value = evaluate(polynomial, pad_to_degree(root), next_precision)
        derivative = polynomial_derivative(polynomial)
        derivative_value = evaluate(derivative, pad_to_degree(root), next_precision)
        derivative_inverse = inverse_unit(derivative_value, next_precision)
        correction = multiply_unit(value, derivative_inverse, next_precision)
        correction = multiply_unit(correction, [modulus - 1], next_precision)
        root = modulo(add(pad_to_degree(root), correction, modulus), polynomial, modulus)
        current_precision = next_precision
      end

      pad_to_degree(root.map { |coefficient| coefficient % prime**target_precision })
    end

    def polynomial_derivative(polynomial_coefficients)
      polynomial_coefficients.each_with_index.each_with_object([]) do |(coefficient, index), derivative|
        derivative << coefficient * index unless index.zero?
      end
    end
  end
end
