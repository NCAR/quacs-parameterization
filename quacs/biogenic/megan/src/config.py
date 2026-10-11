"""Centralized YAML configuration and file path resolution for MEGAN."""
from dataclasses import dataclass
from pathlib import Path
from typing import Any
import yaml

@dataclass(frozen=True)
class MeganConfig:
    root: Path
    data: dict[str, Any]

    @classmethod
    def from_yaml(cls, filename):
        filename = Path(filename).resolve()
        with filename.open(encoding="utf-8") as handle:
            data = yaml.safe_load(handle)
        if not isinstance(data, dict) or not isinstance(data.get("files"), dict) or not isinstance(data.get("model"), dict):
            raise ValueError("YAML requires 'files' and 'model' mappings")
        for name in ("n_class", "n_species", "NLayers", "NRTYP", "NGFORM"):
            value = data["model"].get(name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"model.{name} must be a positive integer")
        mechanism = data["model"].get("mechanism")
        if not isinstance(mechanism, str) or not mechanism.strip():
            raise ValueError("model.mechanism must be a non-empty string (CBMZ or SPECIES)")
        mechanism = mechanism.strip().upper()
        if mechanism not in ("CBMZ", "SPECIES", "NONE"):
            raise ValueError(f"Unsupported model.mechanism={mechanism!r}; supported: CBMZ, SPECIES")
        data["model"]["mechanism"] = mechanism
        sections = data["files"]
        for key in ("inputs", "intermediate_outputs", "outputs"):
            if not isinstance(sections.get(key), dict):
                raise ValueError(f"files.{key} must be a mapping")
        result_file = sections["outputs"].get("result_file")
        if not isinstance(result_file, str) or not result_file.strip():
            raise ValueError("files.outputs.result_file must be a nonempty file path")
        inputs = sections["inputs"]
        mapping = inputs.get("mapping_files")
        for key in ("species_catalog", "class_to_species", "species_to_CBMZ"):
            if not isinstance(mapping, dict) or not isinstance(mapping.get(key), str):
                raise ValueError(f"files.inputs.mapping_files.{key} must be a path")
        return cls(filename.parent, data)

    def path(self, value):
        path = Path(value).expanduser()
        return path if path.is_absolute() else self.root / path

    @property
    def files(self):
        return self.data["files"]

    @property
    def inputs(self):
        return self.files["inputs"]

    @property
    def intermediate_outputs(self):
        return self.files["intermediate_outputs"]

    @property
    def outputs(self):
        return self.files["outputs"]

    @property
    def dimensions(self):
        return self.data["model"]

    @property
    def raw_ef_files(self):
        raw = self.inputs["raw_ef_files"]
        return {"plant_ef": str(self.path(raw["plant_ef"])),
                "speciation": {name: {"file": str(self.path(item["file"])),
                                     "specfrac_column": item["specfrac_column"]}
                               for name, item in raw["speciation"].items()}}

    @property
    def mapping_files(self):
        return {name: self.path(path) for name, path in self.inputs["mapping_files"].items()}
