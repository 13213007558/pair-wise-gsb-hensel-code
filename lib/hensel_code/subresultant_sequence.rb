# frozen_string_literal: true

module HenselCode
  # Subresultant polynomial remainder sequence of two formal polynomials,
  # together with its p-adic precision analysis.
  #
  # The j-th entry is the j-th subresultant polynomial; the 0-th one is the
  # resultant. Each subresultant is known modulo p^precision, where the
  # precision is the smallest precision among the p-adic coefficients of the
  # two operands (infinite when every coefficient is an exact Integer).
  #
  # A subresultant whose coefficients all vanish modulo p^precision is NOT
  # certified to be zero: its coefficients only have a certified valuation
  # lower bound. Over the p-adics, the degree of the gcd of two polynomials is
  # the smallest j with a nonzero j-th subresultant, so the gcd degree can
  # only be certified when the vanishing subresultants are known exactly.
  class SubresultantSequence
    include Enumerable
    include HenselCode::Tools

    attr_reader :prime, :precision

    def initialize(prime:, precision:, integer_polynomials:)
      @prime = prime
      @precision = precision
      @integer_polynomials = integer_polynomials
    end

    def [](index)
      coefficients = @integer_polynomials[index].map { |coefficient| wrap_coefficient(coefficient) }
      FormalPolynomial.new(prime, coefficients)
    end

    def each(&block)
      to_a.each(&block)
    end

    def size
      @integer_polynomials.size
    end

    def to_a
      (0...size).map { |index| self[index] }
    end

    def exact?
      precision.infinite?
    end

    def modulus
      exact? ? nil : prime**precision
    end

    # indices of the subresultants whose coefficients all vanish modulo
    # p^precision; those are undetermined, not certified zeros
    def vanishing_degrees
      (0...size).select { |index| vanishing?(index) }
    end

    # indices of the subresultants certified to be nonzero
    def determined_degrees
      (0...size).reject { |index| vanishing?(index) }
    end

    # certified p-adic valuation of the coefficient of x^i in the j-th
    # subresultant, or nil when it only vanishes to the current precision
    def coefficient_valuation(index, power)
      coefficient = residue(@integer_polynomials[index][power])
      return nil if coefficient.zero? && !exact?

      padic_valuation(coefficient, prime)
    end

    # certified lower bound for the valuation of every coefficient of the
    # j-th subresultant
    def valuation_bound(index)
      nonzero = @integer_polynomials[index].map { |coefficient| residue(coefficient) }.reject(&:zero?)
      return precision_source if nonzero.empty?

      nonzero.map { |coefficient| padic_valuation(coefficient, prime) }.min
    end

    # degree of the p-adic gcd of the two polynomials, or nil when the
    # current precision cannot certify it; a resultant certified nonzero
    # always certifies gcd degree 0, higher degrees require exact vanishing
    # of the leading subresultants
    def gcd_degree
      return 0 unless vanishing?(0)
      return nil unless exact?

      (0...size).find { |index| !vanishing?(index) }
    end

    # certified lower and upper bounds for the p-adic gcd degree
    def gcd_degree_bounds
      lower = exact? ? gcd_degree : 0
      upper = determined_degrees.first || (size - 1)
      [lower, upper]
    end

    private

    def vanishing?(index)
      @integer_polynomials[index].all? { |coefficient| residue(coefficient).zero? }
    end

    def residue(coefficient)
      exact? ? coefficient : coefficient % modulus
    end

    def precision_source
      exact? ? Float::INFINITY : precision
    end

    def wrap_coefficient(coefficient)
      return coefficient if exact?

      # a Rational input keeps the residue intact: constructing a
      # TruncatedFinitePadicExpansion from an Integer triggers a decode/encode
      # roundtrip that can rewrite the residue when the reconstructed
      # denominator is not a p-adic unit
      TruncatedFinitePadicExpansion.new(prime, precision, Rational(coefficient % modulus, 1))
    end
  end
end
