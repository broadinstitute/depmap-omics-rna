from typing import List, Optional

from .base_model import BaseModel


class RnaSequencingAlignments(BaseModel):
    records: List["RnaSequencingAlignmentsRecords"]


class RnaSequencingAlignmentsRecords(BaseModel):
    index_url: Optional[str]
    omics_sequencing_id: str
    reference_genome: Optional[str]
    sequencing_alignment_source: str
    url: str


RnaSequencingAlignments.model_rebuild()
