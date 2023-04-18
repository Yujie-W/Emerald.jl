module Text

using CSV: File
using DataFrames: DataFrame


"""

    read_csv(file::String; skiprows::Int = 0)

Read in the CSV file a data frame, given
- `file` Path to CSV file
- `skiprows` Rows to skip

"""
function read_csv end

read_csv(file::String; skiprows::Int = 0) = DataFrame(File(file; header = skiprows + 1));


end # module
