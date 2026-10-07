# frozen_string_literal: true

module HenselCode
  # Formal polynomial with exact or p-adic coefficients.
  #
  # Unlike HenselCode::Polynomial, which stores the digits of a single p-adic
  # number and is used for carrying in finite p-adic expansions, a
  # FormalPolynomial is an element of a polynomial ring. Each coefficient is
  # either an exact Integer or a p-adic number known to finite precision, i.e.
  # a TruncatedFinitePadicExpansion or a FinitePadicExpansion. Coefficients are
  # given constant term first.
  #
  # The resultant and the subresultants are computed exactly over the integers
  # from representatives of the coefficients and then reduced modulo p^r, where
  # r is the smallest precision among the p-adic coefficients of both operands.
  # This is sound because the resultant and the subresultant coefficients are
  # polynomials with integer coefficients in the coefficients of the inputs.
  #
  # If a leading coefficient vanishes modulo p^r, the degree of the polynomial
  # drops to the current precision, and the computation is carried out with the
  # reduced degree.
  class FormalPolynomial
    attr_reader :prime, :coefficients

    def initialize(prime, coefficients)
      @prime = prime
      @coefficients = coefficients
      valid_prime?
      valid_coefficients?
    end

    def degree
      coefficients.size - 1
    end

    def precision
      exponents = padic_coefficients.map(&:exponent)
      exponents.empty? ? Float::INFINITY : exponents.min
    end

    def exact?
      precision.infinite?
    end

    def resultant(other)
      valid_operand?(other)
      ResultantResult.new(prime: prime, precision: combined_precision(other), value: exact_resultant(other))
    end

    def subresultants(other)
      valid_operand?(other)
      precision = combined_precision(other)
      polys = self.class.integer_subresultants(lifted_coefficients(precision), other.lifted_coefficients(precision))
      SubresultantSequence.new(prime: prime, precision: precision, integer_polynomials: polys)
    end

    # true if the polynomials are known to have a common factor, false if they
    # are known not to have one, and nil if the current precision cannot decide
    def common_factor?(other)
      resultant(other).common_factor?
    end

    # degree of the p-adic gcd, or nil if the current precision cannot decide
    def gcd_degree(other)
      subresultants(other).gcd_degree
    end

    def to_s
      coefficients.map.with_index do |coefficient, index|
        "#{coefficient_string(coefficient)}#{formal_variable(index)}"
      end.join(" + ")
    end

    def inspect
      "<FormalPolynomial: #{self}>"
    end

    # exact resultant and subresultants of integer coefficient arrays
    # (constant term first); used internally and shared with the result classes
    class << self
      def integer_resultant(coefficients1, coefficients2)
        bareiss_determinant(sylvester_matrix(*normalize(coefficients1, coefficients2)))
      end

      def integer_subresultants(coefficients1, coefficients2)
        coefficients1, coefficients2 = normalize(coefficients1, coefficients2)
        swapped = coefficients1.size < coefficients2.size
        coefficients1, coefficients2 = coefficients2, coefficients1 if swapped
        deg1 = coefficients1.size - 1
        deg2 = coefficients2.size - 1
        return [[1]] if deg1.zero? && deg2.zero?

        sequence = (0..deg2).map { |j| subresultant_j(coefficients1, coefficients2, j) }
        return sequence unless swapped

        sequence.map.with_index do |poly, j|
          ((deg1 - j) * (deg2 - j)).odd? ? poly.map { |coefficient| -coefficient } : poly
        end
      end

      private

      def normalize(coefficients1, coefficients2)
        [strip_leading_zeros(coefficients1), strip_leading_zeros(coefficients2)]
      end

      def strip_leading_zeros(coefficients)
        stripped = coefficients.dup
        stripped.pop while stripped.size.positive? && stripped.last.zero?
        raise HenselCode::Error, "polynomial is zero (to the current precision)" if stripped.empty?

        stripped
      end

      def sylvester_matrix(coefficients1, coefficients2)
        deg1 = coefficients1.size - 1
        deg2 = coefficients2.size - 1
        rows = []
        (0...deg2).each { |t| rows << (Array.new(t, 0) + coefficients1.reverse + Array.new(deg2 - 1 - t, 0)) }
        (0...deg1).each { |s| rows << (Array.new(s, 0) + coefficients2.reverse + Array.new(deg1 - 1 - s, 0)) }
        rows
      end

      # j-th subresultant polynomial of f and g with deg f >= deg g, computed
      # as the determinants of the Sylvester submatrices; the 0-th subresultant
      # is exactly the resultant
      def subresultant_j(coefficients1, coefficients2, j)
        deg1 = coefficients1.size - 1
        deg2 = coefficients2.size - 1
        return coefficients2.dup if j == deg1 && j == deg2

        width = deg1 + deg2 - j
        size = deg1 + deg2 - (2 * j)
        rows = []
        (0...(deg2 - j)).each { |t| rows << (Array.new(t, 0) + coefficients1.reverse + Array.new(deg2 - j - 1 - t, 0)) }
        (0...(deg1 - j)).each { |s| rows << (Array.new(s, 0) + coefficients2.reverse + Array.new(deg1 - j - 1 - s, 0)) }
        (0..j).map do |i|
          columns = (0...(size - 1)).to_a << (width - 1 - i)
          bareiss_determinant(rows.map { |row| columns.map { |column| row[column] } })
        end
      end

      # fraction-free Bareiss algorithm: exact determinant of an integer matrix
      def bareiss_determinant(matrix)
        size = matrix.size
        return 1 if size.zero?

        matrix = matrix.map(&:dup)
        sign = 1
        previous = 1
        (0...(size - 1)).each do |k|
          if matrix[k][k].zero?
            pivot = ((k + 1)...size).find { |i| !matrix[i][k].zero? }
            return 0 if pivot.nil?

            matrix[k], matrix[pivot] = matrix[pivot], matrix[k]
            sign = -sign
          end
          ((k + 1)...size).each do |i|
            ((k + 1)...size).each do |j|
              matrix[i][j] = ((matrix[i][j] * matrix[k][k]) - (matrix[i][k] * matrix[k][j])) / previous
            end
          end
          previous = matrix[k][k]
        end
        sign * matrix[size - 1][size - 1]
      end
    end

    protected

    # integer representatives of the coefficients modulo p^precision; exact
    # Integer coefficients are kept as they are
    def lifted_coefficients(precision)
      coefficients.map do |coefficient|
        coefficient.is_a?(Integer) ? coefficient : padic_residue(coefficient) % (prime**precision)
      end
    end

    private

    def valid_prime?
      raise ArgumentError, "prime can't be nil" if prime.nil?
      raise ArgumentError, "prime must be an integer" unless prime.is_a?(Integer)
    end

    def valid_coefficients?
      raise ArgumentError, "coefficients can't be nil" if coefficients.nil?
      raise ArgumentError, "coefficients must be an array" unless coefficients.is_a?(Array)
      raise ArgumentError, "coefficients can't be empty" if coefficients.empty?

      coefficients.each { |coefficient| valid_coefficient?(coefficient) }
    end

    def valid_coefficient?(coefficient)
      return if coefficient.is_a?(Integer)

      unless padic_coefficient?(coefficient)
        raise ArgumentError, "coefficients must be Integers or p-adic Hensel codes"
      end
      if coefficient.prime != prime
        raise HenselCodesWithDifferentPrimes, "coefficient prime differs from polynomial prime"
      end
    end

    def valid_operand?(other)
      raise WrongHenselCodeInputType, "operand must be a FormalPolynomial" unless other.is_a?(FormalPolynomial)
      raise HenselCodesWithDifferentPrimes, "polynomials must have the same prime" if prime != other.prime
    end

    def padic_coefficient?(coefficient)
      !coefficient.is_a?(Integer) && coefficient.respond_to?(:hensel_code) &&
        coefficient.respond_to?(:exponent) && coefficient.respond_to?(:prime)
    end

    def padic_coefficients
      coefficients.select { |coefficient| padic_coefficient?(coefficient) }
    end

    def combined_precision(other)
      [precision, other.precision].min
    end

    def exact_resultant(other)
      precision = combined_precision(other)
      self.class.integer_resultant(lifted_coefficients(precision), other.lifted_coefficients(precision))
    end

    def padic_residue(coefficient)
      digits = coefficient.hensel_code
      return digits if digits.is_a?(Integer)

      digits.each_with_index.sum { |digit, index| digit * (prime**index) }
    end

    def coefficient_string(coefficient)
      coefficient.is_a?(Integer) ? coefficient.to_s : padic_residue(coefficient).to_s
    end

    def formal_variable(index)
      i = index > 1 ? 2 : index
      ["", "x", "x^#{index}"][i]
    end
  end
end
