"""Input readers. v0.2 adds format adapters (MTX, H5AD, 10x H5) here."""

from .readers import (  # noqa: F401
    Matrix,
    ValueState,
    count_lines,
    detect_value_state,
    discover_files,
    guess_sep,
    is_gzip,
    open_text,
    read_dense_matrix,
    read_embedding,
    read_guide_assignments,
    read_header,
    read_table,
)
