# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#               National Center for Biotechnology Information
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government have not placed any restriction on its use or reproduction.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
#  Please cite the author in any work or product based on this material.
#
# ===========================================================================
# Script name: hadoop_json_annotation.py
# Description: a demo script to parse dbSNP RS JSON object and extract the
# records of rs annotation.  The script will produce tab-delimited output
# containing snp_id, seq_id, annotation_release, gene_name, gene_id,
# gene_locus, gene_is_pseudo, gene_orientation, rna_id,
# transcript_change_seq_id, transcript_change_position,
# transcript_change_deleted_sequence, transcript_change_inserted_sequence,
# sequence_ontology_name, sequence_ontology_accession, product_id, spdi_seq_id,
# spdi_position, spdi_deleted_sequence, spdi_inserted_sequence,
# protein_sequence_ontology_name, protein_sequence_ontology_accession.
#
# Sample use:
# python hadoop_json_annotation.py            \
#     -r hadoop hdfs:///path/to/input         \
#     -o hdfs:///path/to/output               \
#     --no-output                             \
#     --jobconf mapreduce.job.name = test.mrjob \
#     --jobconf mapreduce.job.reduces=100     \
#
# Author:  Qiang Wang  wangq2@ncbi.nlm.nih.gov
# For help please contact: tkt-varhd@ncbi.nlm.nih.gov
#
#
# ---------------------------------------------------------------------------

import json

from mrjob.job import MRJob


class MRJsonProcessor(MRJob):

    def mapper(self, _, line):

        data = json.loads(line)

        snp_id = data["refsnp_id"]

        annotations = data["primary_snapshot_data"]["allele_annotations"]

        for a in annotations:
            assembly_annotations = a["assembly_annotation"]
            for aa in assembly_annotations:
                seq_id = aa["seq_id"]
                annotation_release = aa["annotation_release"]
                genes = aa["genes"]

                for g in genes:
                    gene_name = g["name"]
                    gene_id = str(g["id"])
                    gene_locus = g["locus"]
                    gene_is_pseudo = str(int(g["is_pseudo"] == 'true'))
                    gene_orientation = str(int(g["orientation"] == 'minus'))
                    rnas = g["rnas"]

                    for r in rnas:
                        rna_id = str(r["id"])

                        if "transcript_change" in r:
                            transcript_change = r["transcript_change"]
                        else:
                            transcript_change = {"seq_id": "",
                                                 "position": "",
                                                 "deleted_sequence": "",
                                                 "inserted_sequence": ""}

                        sequence_ontology = r["sequence_ontology"]

                        if "product_id" in r:
                            product_id = r["product_id"]
                        else:
                            product_id = ""

                        if "protein" in r:
                            protein = r["protein"]
                            if "spdi" not in protein["variant"]:
                                protein["variant"]["spdi"] = {
                                    "seq_id": "",
                                    "position": "",
                                    "deleted_sequence": "",
                                    "inserted_sequence": ""}
                        else:
                            protein = {"variant":
                                       {"spdi":
                                        {"seq_id": "",
                                         "position": "",
                                         "deleted_sequence": "",
                                         "inserted_sequence": ""
                                         }}, "sequence_ontology": []}

                        prot_spdi = protein["variant"]["spdi"]
                        if len(sequence_ontology) == 0 and \
                                len(protein["sequence_ontology"]) == 0:
                            values = (
                                snp_id, seq_id, annotation_release,
                                gene_name, gene_id, gene_locus,
                                gene_is_pseudo, gene_orientation,
                                rna_id, transcript_change["seq_id"],
                                str(transcript_change["position"]),
                                transcript_change["deleted_sequence"],
                                transcript_change["inserted_sequence"],
                                "", "", product_id,
                                prot_spdi["seq_id"],
                                str(prot_spdi["position"]),
                                prot_spdi["deleted_sequence"],
                                prot_spdi["inserted_sequence"],
                                "", "")
                            print("\t".join(values))
                        elif len(sequence_ontology) != 0 and \
                                len(protein["sequence_ontology"]) == 0:
                            for s in sequence_ontology:
                                values = (
                                    snp_id, seq_id, annotation_release,
                                    gene_name, gene_id, gene_locus,
                                    gene_is_pseudo, gene_orientation,
                                    rna_id, transcript_change["seq_id"],
                                    str(transcript_change["position"]),
                                    transcript_change["deleted_sequence"],
                                    transcript_change["inserted_sequence"],
                                    s["name"], s["accession"], product_id,
                                    prot_spdi["seq_id"],
                                    str(prot_spdi["position"]),
                                    prot_spdi["deleted_sequence"],
                                    prot_spdi["inserted_sequence"],
                                    "", "")
                                print("\t".join(values))
                        elif len(sequence_ontology) == 0 and \
                                len(protein["sequence_ontology"]) != 0:
                            for s in protein["sequence_ontology"]:
                                values = (
                                    snp_id, seq_id, annotation_release,
                                    gene_name, gene_id, gene_locus,
                                    gene_is_pseudo, gene_orientation,
                                    rna_id, transcript_change["seq_id"],
                                    str(transcript_change["position"]),
                                    transcript_change["deleted_sequence"],
                                    transcript_change["inserted_sequence"],
                                    "", "", product_id,
                                    prot_spdi["seq_id"],
                                    str(prot_spdi["position"]),
                                    prot_spdi["deleted_sequence"],
                                    prot_spdi["inserted_sequence"],
                                    s["name"], s["accession"])
                                print("\t".join(values))
                        else:
                            for s in sequence_ontology:
                                for ps in protein["sequence_ontology"]:
                                    values = (
                                        snp_id, seq_id, annotation_release,
                                        gene_name, gene_id, gene_locus,
                                        gene_is_pseudo, gene_orientation,
                                        rna_id, transcript_change["seq_id"],
                                        str(transcript_change["position"]),
                                        transcript_change["deleted_sequence"],
                                        transcript_change["inserted_sequence"],
                                        s["name"], s["accession"], product_id,
                                        prot_spdi["seq_id"],
                                        str(prot_spdi["position"]),
                                        prot_spdi["deleted_sequence"],
                                        prot_spdi["inserted_sequence"],
                                        ps["name"], ps["accession"])
                                    print("\t".join(values))


if __name__ == '__main__':
    MRJsonProcessor.run()
