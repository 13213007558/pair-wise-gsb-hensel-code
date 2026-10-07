# frozen_string_literal: true

module HenselCode
  # Exact integer determinant by fraction-free Bareiss elimination.
  module IntegerDeterminant
    module_function

    def determinant(matrix)
      return 1 if matrix.empty?

      rows = matrix.map(&:dup)
      size = rows.size
      raise ArgumentError, "matrix must be square" unless rows.all? { |row| row.size == size }

      sign = 1
      previous_pivot = 1
      (0...size).each do |pivot_column|
        pivot_row = rows.index.with_index { |_, index| index >= pivot_column && rows[index][pivot_column] != 0 }
        return 0 if pivot_row.nil?

        if pivot_row != pivot_column
          rows[pivot_column], rows[pivot_row] = rows[pivot_row], rows[pivot_column]
          sign = -sign
        end

        pivot = rows[pivot_column][pivot_column]
        (pivot_column + 1...size).each do |row|
          (pivot_column + 1...size).each do |column|
            value = (pivot * rows[row][column]) - (rows[row][pivot_column] * rows[pivot_column][column])
            rows[row][column] = value / previous_pivot
          end
        end
        previous_pivot = pivot
      end

      sign * rows.last.last
    end
  end
end
