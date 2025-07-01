from typing import Any, Dict

from .base_client import BaseClient
from .rna_sequencing_alignments import RnaSequencingAlignments


def gql(q: str) -> str:
    return q


class GumboClient(BaseClient):
    def rna_sequencing_alignments(self, **kwargs: Any) -> RnaSequencingAlignments:
        query = gql(
            """
            query RnaSequencingAlignments {
              records: omics_mapping(where: {datatype: {_eq: "rna"}}) {
                model_id
                model_condition_id
                omics_profile_id
                omics_sequencing_id
                model {
                  cell_line_name
                  stripped_cell_line_name
                }
                omics_sequencing {
                  stranded
                  sequencing_alignments {
                    sequencing_alignment_id: id
                    url
                    index_url
                    reference_genome
                    sequencing_alignment_source
                    size
                  }
                }
              }
            }
            """
        )
        variables: Dict[str, object] = {}
        response = self.execute(
            query=query,
            operation_name="RnaSequencingAlignments",
            variables=variables,
            **kwargs
        )
        data = self.get_data(response)
        return RnaSequencingAlignments.model_validate(data)
