# frozen_string_literal: true

module HenselCode
  # Sparse multivariate polynomial with integer coefficients.
  #
  # Polynomials are given in a structured form: an ordered list of variable
  # names plus a hash mapping exponent arrays to integer coefficients. No
  # string expressions are parsed. For example,
  #
  #   HenselCode::MultivariatePolynomial.new([:x, :y], [2, 0] => 1, [0, 1] => -1)
  #
  # represents the polynomial x^2 - y.
  class MultivariatePolynomial
    attr_reader :variables, :terms

    def initialize(variables, terms)
      @variables = validate_variables(variables)
      @terms = normalize_terms(validate_terms(terms))
    end

    def evaluate(point)
      validate_point(point)
      terms.sum do |exponents, coefficient|
        monomial = variables.each_with_index.map { |name, i| point[name]**exponents[i] }.reduce(1, :*)
        coefficient * monomial
      end
    end

    def partial(variable)
      index = variables.index(variable)
      raise ArgumentError, "unknown variable #{variable.inspect}" if index.nil?

      new_terms = {}
      terms.each do |exponents, coefficient|
        exponent = exponents[index]
        next if exponent.zero?

        new_exponents = exponents.dup
        new_exponents[index] = exponent - 1
        new_terms[new_exponents] = coefficient * exponent
      end
      self.class.new(variables, new_terms)
    end

    def degree
      return 0 if terms.empty?

      terms.keys.map(&:sum).max
    end

    def zero?
      terms.empty?
    end

    def ==(other)
      other.is_a?(self.class) && variables == other.variables && terms == other.terms
    end

    def to_s
      return "0" if terms.empty?

      parts = ordered_terms.map { |exponents, coefficient| term_to_s(exponents, coefficient) }
      first_sign, first_body = parts.first
      rest = parts[1..].map { |sign, body| " #{sign} #{body}" }.join
      "#{first_sign == '-' ? '-' : ''}#{first_body}#{rest}"
    end

    def inspect
      "<MultivariatePolynomial: #{self}>"
    end

    private

    def validate_variables(variables)
      unless variables.is_a?(Array) && !variables.empty? && variables.all? { |v| v.is_a?(Symbol) }
        raise ArgumentError, "variables must be a non-empty array of symbols"
      end
      raise ArgumentError, "variables must be distinct" unless variables.uniq.size == variables.size

      variables.freeze
    end

    def validate_terms(terms)
      raise ArgumentError, "terms must be a hash" unless terms.is_a?(Hash)

      terms.each do |exponents, coefficient|
        valid_exponents = exponents.is_a?(Array) &&
                          exponents.size == variables.size &&
                          exponents.all? { |e| e.is_a?(Integer) && !e.negative? }
        unless valid_exponents
          raise ArgumentError, "term keys must be arrays of #{variables.size} non-negative integer exponents"
        end
        raise ArgumentError, "term coefficients must be integers" unless coefficient.is_a?(Integer)
      end
      terms
    end

    def normalize_terms(terms)
      terms.each_with_object({}) do |(exponents, coefficient), normalized|
        normalized[exponents.freeze] = coefficient unless coefficient.zero?
      end.freeze
    end

    def validate_point(point)
      missing = variables.reject { |name| point.is_a?(Hash) && point[name].is_a?(Integer) }
      unless missing.empty?
        raise ArgumentError, "point must assign an integer to each variable; missing: #{missing.inspect}"
      end
    end

    def ordered_terms
      terms.sort_by { |exponents, _| [-exponents.sum, exponents.map { |e| -e }] }
    end

    def term_to_s(exponents, coefficient)
      monomial = monomial_to_s(exponents)
      return [coefficient.negative? ? "-" : "+", coefficient.abs.to_s] if monomial.empty?

      coefficient_string = coefficient.abs == 1 ? "" : "#{coefficient.abs}*"
      [coefficient.negative? ? "-" : "+", "#{coefficient_string}#{monomial}"]
    end

    def monomial_to_s(exponents)
      variables.zip(exponents).select { |_, exponent| exponent.positive? }.map do |name, exponent|
        exponent == 1 ? name.to_s : "#{name}^#{exponent}"
      end.join("*")
    end
  end
end
