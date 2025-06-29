"""
Gene ID renaming module for EVM post-processing.
Renames EVM-generated gene IDs to user-defined prefix format.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
from collections import defaultdict
import logging

from .species_config import GeneNamingConfig


logger = logging.getLogger(__name__)


class GeneIDRenamer:
    """
    Renames EVM-generated gene IDs to user-defined prefix format.
    
    Converts EVM IDs like 'evm.model.Chr5.1' to prefix-based format
    like 'NBS_0001' or 'CaCM334_0001'.
    """
    
    def __init__(self, prefix: str = None):
        """
        Initialize the gene ID renamer.
        
        Args:
            prefix: User-defined prefix for gene IDs (default: "NBS")
        """
        self.naming_config = GeneNamingConfig(prefix)
        self.gene_counter = 0  # Global gene counter
        self.feature_counters = {}  # Track feature counters per gene
        self.id_mapping = {}  # Old ID -> New ID mapping
        self.reverse_mapping = {}  # New ID -> Old ID mapping
        self.gene_associations = {}  # Track gene-transcript associations
        
        logger.info(f"Initialized GeneIDRenamer with prefix: {self.naming_config.get_prefix()}")
    
    def _parse_evm_id(self, evm_id: str) -> Optional[Tuple[str, str, str, str]]:
        """
        Parse EVM ID to extract components.
        
        Expected formats:
        - evm.model.Chr5.1
        - evm.TU.Chr5.1
        - cds.evm.model.Chr5.1.1 (CDS with additional number)
        
        Args:
            evm_id: Original EVM ID
            
        Returns:
            Tuple of (prefix, feature_type, chromosome, number) or None
        """
        patterns = [
            r'^cds\.evm\.(model|TU)\.(.+?)\.(\d+)\.(\d+)$',     # CDS with suffix: cds.evm.model.Chr.N.X
            r'^cds\.evm\.(model|TU)\.(.+?)\.(\d+)$',            # CDS without suffix: cds.evm.model.Chr.N
            r'^evm\.(model|TU)\.(.+?)\.(\d+)$',                 # Standard mRNA format: evm.model.Chr.N
        ]
        
        for pattern in patterns:
            match = re.match(pattern, evm_id)
            if match:
                groups = match.groups()
                if len(groups) == 4:  # CDS format with suffix: cds.evm.model.Chr.N.X
                    feature_type, chromosome, number, cds_suffix = groups
                    return ("cds.", feature_type, chromosome, number)
                elif len(groups) == 3:  # CDS format without suffix OR standard mRNA format
                    # Check if this is a CDS (starts with cds.) or mRNA (starts with evm.)
                    if evm_id.startswith('cds.'):
                        feature_type, chromosome, number = groups
                        return ("cds.", feature_type, chromosome, number)
                    else:
                        feature_type, chromosome, number = groups
                        return ("", feature_type, chromosome, number)
        
        return None
    
    def _get_next_gene_number(self) -> int:
        """
        Get next sequential gene number.
        
        Returns:
            Next gene number
        """
        self.gene_counter += 1
        return self.gene_counter
    
    def _get_base_evm_id(self, evm_id: str) -> str:
        """
        Get base EVM ID for grouping related features.
        
        For CDS features like 'cds.evm.model.Chr1.1.1', this returns
        'evm.model.Chr1.1' to group them with their parent mRNA.
        For exon features like 'evm.model.Chr1.1.exon1', this also returns
        'evm.model.Chr1.1' to group them with their parent mRNA.
        
        Args:
            evm_id: Original EVM ID
            
        Returns:
            Base EVM ID for grouping
        """
        # Handle exon features with pattern: evm.model.Chr1.1.exon1
        if '.exon' in evm_id:
            # Extract the base part before .exon
            base_part = evm_id.split('.exon')[0]
            return base_part
        
        parsed = self._parse_evm_id(evm_id)
        if parsed:
            prefix, feature_type, chromosome, number = parsed
            # For CDS features, strip the final .X suffix to group with parent mRNA
            if prefix and prefix.startswith("cds"):
                # cds.evm.model.Chr1.1.1 -> evm.model.Chr1.1
                return f"evm.model.{chromosome}.{number}"
            else:
                # evm.model.Chr1.1 -> evm.model.Chr1.1
                return f"evm.model.{chromosome}.{number}"
        return evm_id
    
    def _get_next_feature_number(self, base_evm_id: str, feature_type: str) -> int:
        """
        Get next feature number for a specific gene and feature type.
        
        Args:
            base_evm_id: Base EVM ID for grouping
            feature_type: Feature type (exon, CDS, etc.)
            
        Returns:
            Next feature number
        """
        key = f"{base_evm_id}_{feature_type}"
        if key not in self.feature_counters:
            self.feature_counters[key] = 0
        self.feature_counters[key] += 1
        return self.feature_counters[key]
    
    def _generate_new_id(self, old_id: str, feature_type: str) -> str:
        """
        Generate new prefix-based gene ID with appropriate suffix.
        
        Args:
            old_id: Original EVM ID
            feature_type: GFF3 feature type (gene, mRNA, exon, CDS, etc.)
            
        Returns:
            New prefix-based gene ID
        """
        # 获取基础EVM ID用于分组
        base_evm_id = self._get_base_evm_id(old_id)
        
        # 如果这个基础ID还没有分配基因编号，分配一个
        if base_evm_id not in self.gene_associations:
            gene_number = self._get_next_gene_number()
            self.gene_associations[base_evm_id] = gene_number
        else:
            gene_number = self.gene_associations[base_evm_id]
        
        # 根据特征类型生成相应的ID
        if feature_type.lower() in ['gene', 'mrna']:
            new_id = self.naming_config.generate_gene_id(gene_number, feature_type)
        else:
            # 对于其他特征类型（如exon, CDS），使用基因编号加特征编号
            feature_number = self._get_next_feature_number(base_evm_id, feature_type.lower())
            prefix = self.naming_config.get_prefix()
            new_id = f"{prefix}_{gene_number:04d}.{feature_type.lower()}{feature_number}"
        
        logger.debug(f"Mapped {old_id} -> {new_id} (type: {feature_type}, base: {base_evm_id})")
        return new_id
    
    def rename_gff3_file(self, input_file: str, output_file: str) -> Dict[str, str]:
        """
        Rename gene IDs in GFF3 file.
        
        Args:
            input_file: Input GFF3 file with EVM IDs
            output_file: Output GFF3 file with renamed IDs
            
        Returns:
            Dictionary mapping old IDs to new IDs
        """
        logger.info(f"Renaming gene IDs in {input_file} -> {output_file}")
        
        self.id_mapping.clear()
        self.reverse_mapping.clear()
        self.gene_counter = 0
        self.gene_associations.clear()
        self.feature_counters.clear()
        self._used_cds_keys = set()  # Reset CDS key tracking
        
        # Store mRNA information for gene creation
        mrna_features = []
        
        # First pass: collect all EVM IDs and store mRNA features
        with open(input_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                seqid, source, feature_type = fields[0], fields[1], fields[2]
                start, end, score, strand, phase = fields[3], fields[4], fields[5], fields[6], fields[7]
                attributes = fields[8]
                
                # Extract ID from attributes
                id_match = re.search(r'ID=([^;]+)', attributes)
                if id_match:
                    old_id = id_match.group(1)
                    
                    # Only process EVM IDs
                    parsed = self._parse_evm_id(old_id)
                    if parsed:
                        # For CDS features with duplicate IDs, create unique mappings
                        if feature_type.lower() == 'cds' and old_id.startswith('cds.'):
                            # Create a unique key for each CDS occurrence
                            base_evm_id = self._get_base_evm_id(old_id)
                            unique_key = f"{old_id}#{len([k for k in self.id_mapping.keys() if k.startswith(old_id)])}"
                            if unique_key not in self.id_mapping:
                                new_id = self._generate_new_id(old_id, feature_type)
                                self.id_mapping[unique_key] = new_id
                                self.reverse_mapping[new_id] = unique_key
                        else:
                            # Generate new ID if not already mapped
                            if old_id not in self.id_mapping:
                                new_id = self._generate_new_id(old_id, feature_type)
                                self.id_mapping[old_id] = new_id
                                self.reverse_mapping[new_id] = old_id
                        
                        # Store mRNA features for later gene creation
                        if feature_type.lower() == 'mrna':
                            mrna_features.append({
                                'seqid': seqid,
                                'source': source,
                                'start': int(start),
                                'end': int(end),
                                'score': score,
                                'strand': strand,
                                'phase': phase,
                                'old_id': old_id,
                                'new_id': self.id_mapping[old_id],
                                'base_evm_id': self._get_base_evm_id(old_id)
                            })
        
        # Second pass: organize and write output by gene groups
        self._write_output_by_gene_groups(input_file, output_file, mrna_features)
        
        logger.info(f"Renamed {len(self.id_mapping)} gene IDs")
        return self.id_mapping.copy()
    
    def _write_output_by_gene_groups(self, input_file: str, output_file: str, mrna_features: List[Dict]) -> None:
        """
        Write output file organized by gene groups.
        Each gene and all its related features (mRNA, exon, CDS) are grouped together.
        
        Args:
            input_file: Input GFF3 file path
            output_file: Output GFF3 file path
            mrna_features: List of mRNA feature information
        """
        # Read all lines from input file and organize by base EVM ID
        all_lines = []
        with open(input_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                all_lines.append(line.strip())
        
        # Group features by base EVM ID
        features_by_gene = defaultdict(list)
        
        for line in all_lines:
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            
            feature_type = fields[2]
            attributes = fields[8]
            
            # Skip existing gene features
            if feature_type.lower() == 'gene':
                continue
            
            # Extract ID and determine base EVM ID for grouping
            id_match = re.search(r'ID=([^;]+)', attributes)
            if id_match:
                feature_id = id_match.group(1)
                base_evm_id = self._get_base_evm_id(feature_id)
                features_by_gene[base_evm_id].append(line)
        
        # Sort gene groups by gene number for consistent output
        sorted_genes = sorted(features_by_gene.keys(), 
                            key=lambda x: self.gene_associations.get(x, 999))
        
        # Write output file
        with open(output_file, 'w') as outfile:
            # Write GFF3 header
            outfile.write("##gff-version 3\n")
            
            # Process each gene group
            for base_evm_id in sorted_genes:
                # Create and write gene feature
                gene_number = self.gene_associations[base_evm_id]
                gene_id = self.naming_config.generate_gene_id(gene_number, 'gene')
                
                # Find mRNA info for this gene to get coordinates
                gene_mrnas = [m for m in mrna_features if m['base_evm_id'] == base_evm_id]
                if gene_mrnas:
                    mrna = gene_mrnas[0]  # Use first mRNA for gene template
                    gene_start = min(m['start'] for m in gene_mrnas)
                    gene_end = max(m['end'] for m in gene_mrnas)
                    
                    # Write gene feature
                    gene_line = f"{mrna['seqid']}\t{mrna['source']}\tgene\t{gene_start}\t{gene_end}\t{mrna['score']}\t{mrna['strand']}\t.\tID={gene_id}\n"
                    outfile.write(gene_line)
                
                # Write all related features for this gene
                for line in features_by_gene[base_evm_id]:
                    modified_line = self._rename_line_ids(line)
                    outfile.write(modified_line)
                
                # Add blank line between gene groups for readability
                outfile.write("\n")
    
    def _rename_line_ids(self, line: str) -> str:
        """
        Rename IDs in a single GFF3 line.
        
        Args:
            line: GFF3 line
            
        Returns:
            Line with renamed IDs
        """
        fields = line.strip().split('\t')
        if len(fields) < 9:
            return line
        
        feature_type = fields[2]
        attributes = fields[8]
        
        # Replace ID attribute
        def replace_id(match):
            old_id = match.group(1)
            # For CDS features, we need to find the right unique key
            if feature_type.lower() == 'cds' and old_id.startswith('cds.'):
                # Find the next unused mapping for this CDS ID
                for key, new_id in self.id_mapping.items():
                    if key.startswith(f"{old_id}#") and key not in getattr(self, '_used_cds_keys', set()):
                        if not hasattr(self, '_used_cds_keys'):
                            self._used_cds_keys = set()
                        self._used_cds_keys.add(key)
                        return f"ID={new_id}"
                # Fallback to original ID if no mapping found
                return f"ID={old_id}"
            else:
                new_id = self.id_mapping.get(old_id, old_id)
                return f"ID={new_id}"
        
        attributes = re.sub(r'ID=([^;]+)', replace_id, attributes)
        
        # Replace Parent attribute with special handling for mRNA
        def replace_parent(match):
            old_parent = match.group(1)
            
            # For mRNA features, Parent should point to the gene
            if feature_type.lower() == 'mrna':
                # Find the corresponding gene ID
                parsed = self._parse_evm_id(old_parent)
                if parsed:
                    base_evm_id = self._get_base_evm_id(old_parent)
                    if base_evm_id in self.gene_associations:
                        gene_number = self.gene_associations[base_evm_id]
                        gene_id = self.naming_config.generate_gene_id(gene_number, 'gene')
                        return f"Parent={gene_id}"
            
            # For other features, use normal mapping
            new_parent = self.id_mapping.get(old_parent, old_parent)
            return f"Parent={new_parent}"
        
        attributes = re.sub(r'Parent=([^;]+)', replace_parent, attributes)
        
        # Replace Name attribute if it matches an old ID
        def replace_name(match):
            old_name = match.group(1)
            new_name = self.id_mapping.get(old_name, old_name)
            return f"Name={new_name}"
        
        attributes = re.sub(r'Name=([^;]+)', replace_name, attributes)
        
        fields[8] = attributes
        return '\t'.join(fields) + '\n'
    
    def get_mapping_summary(self) -> Dict[str, int]:
        """
        Get summary of ID mappings.
        
        Returns:
            Dictionary with mapping statistics
        """
        return {
            "total_genes_renamed": len(self.id_mapping),
            "prefix_used": self.naming_config.get_prefix()
        }
    
    def save_mapping_file(self, output_file: str) -> None:
        """
        Save ID mapping to file for reference.
        
        Args:
            output_file: Output file path
        """
        with open(output_file, 'w') as f:
            f.write("# Gene ID Mapping File\n")
            f.write("# Original_ID\tNew_ID\n")
            for old_id, new_id in sorted(self.id_mapping.items()):
                f.write(f"{old_id}\t{new_id}\n")
        
        logger.info(f"ID mapping saved to {output_file}")
    
    def get_chromosome_summary(self, gff3_file: str) -> Dict[str, int]:
        """
        Get summary of chromosomes in GFF3 file.
        
        Args:
            gff3_file: GFF3 file to analyze
            
        Returns:
            Dictionary with chromosome names and feature counts
        """
        chromosome_counts = defaultdict(int)
        
        with open(gff3_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) >= 1:
                    seqid = fields[0]
                    chromosome_counts[seqid] += 1
        
        return dict(chromosome_counts)
    
    @classmethod
    def rename_evm_output(cls, input_file: str, output_file: str, 
                         prefix: str = None) -> Dict[str, str]:
        """
        Convenience method to rename EVM output in one call.
        
        Args:
            input_file: Input GFF3 file with EVM IDs
            output_file: Output GFF3 file with renamed IDs
            prefix: User-defined prefix for gene IDs
            
        Returns:
            Dictionary mapping old IDs to new IDs
        """
        renamer = cls(prefix)
        return renamer.rename_gff3_file(input_file, output_file)