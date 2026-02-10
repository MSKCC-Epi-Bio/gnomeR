# Check a Data Frame for Required Columns

Check a Data Frame for Required Columns

## Usage

``` r
.check_required_cols(data, required_cols, add_to_message = NULL)
```

## Arguments

- data:

  A data frame to check

- required_cols:

  A character specifying names of columns to check

- add_to_message:

  a vector (preferrably named) of text to add to the error message for
  specific cases

## Value

If data set doesn't have required columns it will return an error
message. If it does have required columns, nothing will be returned
