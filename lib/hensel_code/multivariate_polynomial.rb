# frozen_string_literal: true

module HenselCode
  # Structured multivariate polynomial with integer coefficients.
  #
  # A polynomial is described by its +arity+ (number of variables) and a
  # sparse +terms+ hash mapping exponent vectors to integer coefficients:
  #
  #   # 3x^2y - 2y + 7 in two variables
  #   HenselCode::MultivariatePolynomial.new(2, [2, 1] => 3, [0, 1] => -2, [0, 0] => 7)
  #
  # Equations are given directly as these objects; no string expressions
  # are parsed or evaluated.
  class MultivariatePolynomial
    attr_reader :arity, :terms

    def initialize(arity, terms)
      @arity = arity
      @terms = terms || {}
      valid_arity?
      valid_terms?
      @terms = normalize_terms(@terms)
    end

    def evaluate(point)
      valid_point?(point)
      terms.sum do |exponents, coefficient|
        coefficient * exponents.zip(point).reduce(1) { |acc, (e, x)| acc * (x**e) }
      end
    end

    def evaluate_mod(point, modulus)
      evaluate(point) % modulus
    end

    def partial(variable_index)
      raise ArgumentError, "variable index out of range" unless variable_index.between?(0, arity - 1)

      new_terms = {}
      terms.each do |exponents, coefficient|
        exponent = exponents[variable_index]
        next if exponent.zero?

        new_exponents = exponents.dup
        new_exponents[variable_index] -= 1
        new_terms[new_exponents] = coefficient * exponent
      end
      self.class.new(arity, new_terms)
    end

    def zero?
      terms.empty?
    end

    def to_s
      return "0" if zero?

      terms.map do |exponents, coefficient|
        variables = exponents.each_with_index.map do |exponent, index|
          next if exponent.zero?

          exponent == 1 ? "x#{index}" : "x#{index}^#{exponent}"
        end.compact
        ([coefficient] + variables).join("*")
      end.join(" + ")
    end

    def inspect
      "<MultivariatePolynomial: #{self}>"
    end

    private

    def normalize_terms(terms)
      terms.each_with_object({}) do |(exponents, coefficient), hash|
        hash[exponents.dup] = coefficient unless coefficient.zero?
      end
    end

    def valid_arity?
      raise ArgumentError, "arity can't be nil" if arity.nil?
      raise ArgumentError, "arity must be an integer" unless arity.is_a?(Integer)
      raise ArgumentError, "arity must be positive" unless arity.positive?
    end

    def valid_terms?
      terms.each do |exponents, coefficient|
        unless exponents.is_a?(Array) && exponents.size == arity && exponents.all? { |e| e.is_a?(Integer) && e >= 0 }
          raise ArgumentError, "term exponents must be arrays of #{arity} non-negative integers"
        end
        raise ArgumentError, "term coefficients must be integers" unless coefficient.is_a?(Integer)
      end
    end

    def valid_point?(point)
      return if point.is_a?(Array) && point.size == arity && point.all? { |x| x.is_a?(Integer) }

      raise ArgumentError, "point must be an array of #{arity} integers"
    end
  end
end
