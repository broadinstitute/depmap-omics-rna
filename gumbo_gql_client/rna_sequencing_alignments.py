from typing import List, Optional

from .base_model import BaseModel


class RnaSequencingAlignments(BaseModel):
    records: List["RnaSequencingAlignmentsRecords"]


class RnaSequencingAlignmentsRecords(BaseModel):
    model_id: Optional[str]
    model_condition_id: Optional[str]
    omics_profile_id: Optional[str]
    omics_sequencing_id: Optional[str]
    model: Optional["RnaSequencingAlignmentsRecordsModel"]
    omics_sequencing: Optional["RnaSequencingAlignmentsRecordsOmicsSequencing"]


class RnaSequencingAlignmentsRecordsModel(BaseModel):
    cell_line_name: Optional[str]
    stripped_cell_line_name: Optional[str]


class RnaSequencingAlignmentsRecordsOmicsSequencing(BaseModel):
    stranded: bool
    sequencing_alignments: List[
        "RnaSequencingAlignmentsRecordsOmicsSequencingSequencingAlignments"
    ]


class RnaSequencingAlignmentsRecordsOmicsSequencingSequencingAlignments(BaseModel):
    sequencing_alignment_id: int
    url: str
    index_url: Optional[str]
    reference_genome: Optional[str]
    sequencing_alignment_source: str
    size: int


RnaSequencingAlignments.model_rebuild()
RnaSequencingAlignmentsRecords.model_rebuild()
RnaSequencingAlignmentsRecordsOmicsSequencing.model_rebuild()
