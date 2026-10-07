# frozen_string_literal: true

module HenselCode
  # Resultants and subresultants for formal polynomials over Z/p^nZ.
  class PolynomialResultant
    include IntegerDeterminant

    attr_reader :f, :g

    def initialize(f, g)
      unless f.is_a?(FormalPolynomial) && g.is_a?(FormalPolynomial)
        raise ArgumentError, "arguments must be FormalPolynomial objects"
      end
      raise ArgumentError, "polynomials must have the same prime and precision" unless compatible?(f, g)

      @f = f
      @g = g
    end

    def self.resultant(f, g)
      new(f, g).resultant
    end

    def self.subresultants(f, g)
      new(f, g).subresultants
    end

    def self.compute(f, g)
      new(f, g).compute
    end

    def resultant
      return scalar_result(0) if f.coefficients.empty? || g.coefficients.empty?

      scalar_result(determinant(sylvester_matrix))
    end

    def subresultants
      return [] if f.coefficients.empty? || g.coefficients.empty?

      (0...[effective_degree(f), effective_degree(g)].min).map do |degree|
        subresultant(degree)
      end
    end

    def subresultant(degree)
      max_degree = [effective_degree(f), effective_degree(g)].min
      raise ArgumentError, "degree must be between 0 and #{max_degree - 1}" if degree.negative? || degree >= max_degree

      matrix_rows = subresultant_rows(degree)
      width = matrix_rows.first.size
      coefficient_columns = ((degree + 1)...width).to_a
      coefficients = (0..degree).map do |power|
        columns = coefficient_columns + [power]
        determinant(matrix_rows.map { |row| columns.map { |column| row[column] } })
      end
      polynomial_result(FormalPolynomial.new(f.prime, f.precision, coefficients))
    end

    def compute
      resultant_report = resultant
      {
        resultant: resultant_report,
        subresultants: subresultants
      }
    end

    private

    def compatible?(left, right)
      left.prime == right.prime && left.precision == right.precision
    end

    def effective_degree(polynomial)
      polynomial.coefficients.empty? ? -1 : polynomial.degree_bound
    end

    def sylvester_matrix
      subresultant_rows(0)
    end

    def subresultant_rows(degree)
      m = effective_degree(f)
      n = effective_degree(g)
      width = m + n - degree
      rows = []
      append_shifted_rows(rows, f, n - degree, width)
      append_shifted_rows(rows, g, m - degree, width)
      rows
    end

    def append_shifted_rows(rows, polynomial, count, width)
      count.times do |shift|
        row = Array.new(width, 0)
        polynomial.coefficients.each_with_index { |coefficient, degree| row[degree + shift] = coefficient }
        rows << row
      end
    end

    def scalar_result(value)
      residue = PadicResidue.new(f.prime, f.precision, value)
      {
        prime: residue.prime,
        precision: residue.precision,
        residue: residue.residue,
        valuation: residue.valuation_determined? ? residue.valuation : nil,
        valuation_lower_bound: residue.valuation,
        valuation_determined: residue.valuation_determined?,
        known_unit_precision: residue.known_precision,
        unresolved_precision: residue.unresolved_precision
      }
    end

    def polynomial_result(polynomial)
      {
        degree_bound: polynomial.degree_bound,
        degree: polynomial.degree,
        coefficients: polynomial.valuation_report
      }
    end
  end
end
