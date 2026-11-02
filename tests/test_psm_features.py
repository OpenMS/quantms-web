import pyopenms as poms

from src.common.psm_features import features_by_id, features_frame, hit_features


def _hit(meta: dict) -> poms.PeptideHit:
    hit = poms.PeptideHit()
    hit.setSequence(poms.AASequence.fromString("PEPTIDEK"))
    hit.setCharge(2)
    hit.setScore(0.001)
    for key, value in meta.items():
        hit.setMetaValue(key, value)
    return hit


COMET_AND_PERCOLATOR = {
    "COMET:xcorr": 2.5, "MS:1002252": 2.5, "MS:1002257": 0.001, "expect": 0.001,
    "MS:1001492": 1.2, "MS:1001491": 0.004, "MS:1001493": 0.01,
    "target_decoy": "target", "protein_references": "unique", "my_new_feature": 3.0,
}


def test_hit_features_names_dedupes_and_orders():
    rows = hit_features(_hit(COMET_AND_PERCOLATOR))
    labels = [r["feature"] for r in rows]
    assert labels[:3] == ["Percolator score", "q-value", "Posterior error probability (PEP)"]
    assert labels.count("XCorr") == 1  # COMET:xcorr is the same value as MS:1002252
    assert "E-value" in labels and "expect" not in labels
    assert "protein_references" not in labels
    # an unknown feature keeps its key, last
    assert rows[-1] == {"group": "Other", "feature": "my_new_feature", "value": 3.0, "help": ""}


def test_features_by_id_follows_parse_idxml_numbering(tmp_path):
    peptides = []
    for i, xcorr in enumerate((1.0, 2.0, 3.0)):
        pep = poms.PeptideIdentification()
        pep.setIdentifier("S")
        pep.setScoreType("expect")
        pep.setHits([_hit({"MS:1002252": xcorr})])
        peptides.append(pep)
    prot = poms.ProteinIdentification()
    prot.setIdentifier("S")
    if hasattr(poms, "PeptideIdentificationList"):
        container = poms.PeptideIdentificationList()
        for pep in peptides:
            container.push_back(pep)
        peptides = container
    path = tmp_path / "run.idXML"
    poms.IdXMLFile().store(str(path), [prot], peptides)

    features = features_by_id(path)
    assert sorted(features) == [0, 1, 2]
    assert [f[0]["value"] for f in (features[0], features[1], features[2])] == [1.0, 2.0, 3.0]


def test_features_frame_shows_values_as_text():
    frame = features_frame(hit_features(_hit({"MS:1001491": 0.0001234, "MS:1002258": 9})))
    assert list(frame.columns) == ["Source", "Feature", "Value", "Meaning"]
    assert frame["Value"].tolist() == ["0.0001234", "9"]
