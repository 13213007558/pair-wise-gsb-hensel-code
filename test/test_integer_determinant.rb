# frozen_string_literal: true

require "test_helper"

class TestIntegerDeterminant < Minitest::Test
  def test_determinant_with_row_exchange
    matrix = [
      [0, 1, 2],
      [3, 0, 1],
      [2, 1, 0]
    ]

    assert_equal 8, HenselCode::IntegerDeterminant.determinant(matrix)
  end

  def test_singular_determinant
    matrix = [
      [1, 2, 3],
      [4, 5, 6],
      [7, 8, 9]
    ]

    assert_equal 0, HenselCode::IntegerDeterminant.determinant(matrix)
  end
end
