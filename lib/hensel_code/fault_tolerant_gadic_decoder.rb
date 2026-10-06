# frozen_string_literal: true

module HenselCode
  # Report produced by HenselCode::FaultTolerantGadicDecoder.
  #
  # +status+ is one of:
  # * +:ok+ - +rational+ is the unique rational within the given bounds that can
  #   be certified from the input. +inconsistent_moduli+ lists the moduli whose
  #   residues disagree with it.
  # * +:ambiguous+ - more than one rational within the bounds explains the
  #   residues, and none of them satisfies the certification condition.
  #   +candidates+ lists all of them.
  # * +:insufficient_information+ - no rational within the bounds could be
  #   certified (either none was reconstructed or the single reconstructed
  #   candidate fails the certification condition).
  #
  # The certificate is independently checkable: the decode is provably unique
  # whenever +consistent_product > uniqueness_bound+, where +consistent_product+
  # is the product of the moduli consistent with +rational+ and
  # +uniqueness_bound+ is +2 * numerator_bound * denominator_bound *
  # error_product+.
  FaultTolerantDecodeResult = Struct.new(
    :status, :rational, :consistent_moduli, :inconsistent_moduli,
    :consistent_product, :error_product, :uniqueness_bound, :candidates,
    keyword_init: true
  ) do
    def ok?
      status == :ok
    end

    def verified?
      !consistent_product.nil? && consistent_product > uniqueness_bound
    end
  end

  # Fault-tolerant decoding of g-adic (CRT-based) Hensel codes.
  #
  # In a distributed computation each residue of a truncated g-adic Hensel code
  # is computed independently, so a few residues may come back corrupted. Given
  # pairwise coprime moduli, the observed residues, and bounds N (numerator) and
  # D (denominator) for the encoded rational, this decoder enumerates the
  # subsets of moduli obtained by dropping up to +max_errors+ of them, combines
  # each kept subset with CRT, and runs rational reconstruction on the combined
  # residue. Each reconstructed rational a/b is then checked against *every*
  # modulus, producing its consistency set S = { i : b * residue_i == a (mod
  # m_i) }.
  #
  # A candidate is certified (provably the unique rational within the bounds)
  # iff
  #
  #   product of moduli in S  >  2 * N * D * product of moduli outside S
  #
  # Proof sketch: if a/b and a'/b' both lie within the bounds and both satisfy
  # the condition with consistency sets S and S', then the product over S and
  # S' jointly exceeds 2ND, so a * b' == a' * b (mod that product) with
  # |a * b' - a' * b| <= 2ND forces a * b' == a' * b, i.e. a/b == a'/b'. Hence
  # no majority voting is involved and no candidate is accepted merely for
  # being small: acceptance requires the checkable product condition above.
  module FaultTolerantGadicDecoder
    extend Tools

    MAX_ERRORS = 2

    def self.decode(moduli, residues, numerator_bound, denominator_bound, max_errors: MAX_ERRORS)
      validate_inputs!(moduli, residues, numerator_bound, denominator_bound)
      total_product = moduli.inject(:*)
      threshold_factor = 2 * numerator_bound * denominator_bound
      candidates = {}

      kept_subsets(moduli.size, max_errors).each do |kept|
        kept_moduli = kept.map { |i| moduli[i] }
        kept_product = kept_moduli.inject(:*)
        next if kept_product <= threshold_factor

        reconstruction = rational_reconstruct(
          crt(kept_moduli, kept.map { |i| residues[i] }), kept_product, numerator_bound, denominator_bound
        )
        next if reconstruction.nil?

        rational = Rational(*reconstruction)
        next if candidates.key?(rational)

        candidates[rational] = candidate_info(moduli, residues, reconstruction, total_product, threshold_factor)
      end

      build_result(moduli, candidates, threshold_factor)
    end

    def self.rational_reconstruct(value, modulus, numerator_bound, denominator_bound)
      _, _, _, numerator, denominator = eea_core(modulus, value % modulus, numerator_bound)
      return nil if denominator.zero?

      numerator = -numerator if denominator.negative?
      denominator = -denominator if denominator.negative?
      return nil if denominator > denominator_bound
      return nil unless numerator.gcd(denominator) == 1 && denominator.gcd(modulus) == 1

      [numerator, denominator]
    end
    private_class_method :rational_reconstruct

    def self.candidate_info(moduli, residues, reconstruction, total_product, threshold_factor)
      numerator, denominator = reconstruction
      consistent = moduli.each_index.select do |i|
        ((denominator * residues[i]) - numerator) % moduli[i] == 0
      end
      consistent_product = consistent.map { |i| moduli[i] }.inject(:*)
      error_product = total_product / consistent_product
      {
        consistent: consistent,
        consistent_product: consistent_product,
        error_product: error_product,
        provable: consistent_product > threshold_factor * error_product
      }
    end
    private_class_method :candidate_info

    def self.build_result(moduli, candidates, threshold_factor)
      provable = candidates.select { |_, info| info[:provable] }

      if provable.keys.size == 1
        ok_result(moduli, provable, candidates, threshold_factor)
      else
        status = candidates.size > 1 ? :ambiguous : :insufficient_information
        FaultTolerantDecodeResult.new(status: status, candidates: candidates.keys)
      end
    end
    private_class_method :build_result

    def self.ok_result(moduli, provable, candidates, threshold_factor)
      rational, info = provable.first
      inconsistent_indexes = moduli.each_index.to_a - info[:consistent]
      FaultTolerantDecodeResult.new(
        status: :ok,
        rational: rational,
        consistent_moduli: info[:consistent].map { |i| moduli[i] },
        inconsistent_moduli: inconsistent_indexes.map { |i| moduli[i] },
        consistent_product: info[:consistent_product],
        error_product: info[:error_product],
        uniqueness_bound: threshold_factor * info[:error_product],
        candidates: candidates.keys
      )
    end
    private_class_method :ok_result

    def self.kept_subsets(size, max_errors)
      indexes = (0...size).to_a
      (0..[max_errors, size - 1].min).flat_map { |drops| indexes.combination(size - drops).to_a }
    end
    private_class_method :kept_subsets

    def self.validate_inputs!(moduli, residues, numerator_bound, denominator_bound)
      unless moduli.is_a?(Array) && residues.is_a?(Array) && !moduli.empty? && moduli.size == residues.size
        raise ArgumentError, "moduli and residues must be non-empty arrays of the same size"
      end
      unless moduli.all? { |modulus| modulus.is_a?(Integer) && modulus > 1 }
        raise ArgumentError, "moduli must be integers greater than 1"
      end
      valid_residues = residues.each_with_index.all? do |residue, i|
        residue.is_a?(Integer) && residue >= 0 && residue < moduli[i]
      end
      unless valid_residues
        raise ArgumentError, "each residue must be an integer in the range of its modulus"
      end
      unless numerator_bound.is_a?(Integer) && numerator_bound >= 1 &&
             denominator_bound.is_a?(Integer) && denominator_bound >= 1
        raise ArgumentError, "numerator and denominator bounds must be positive integers"
      end
      return if moduli.combination(2).all? { |mod1, mod2| mod1.gcd(mod2) == 1 }

      raise NonCoprimeModuli, "moduli must be pairwise coprime"
    end
    private_class_method :validate_inputs!
  end
end
