import os
import pytest
from unittest.mock import patch, MagicMock
from cnvpytor.fasta import Fasta

TEST_FASTA = os.path.join(os.path.dirname(__file__), '../test_data/test.fa.gz')

def test_get_chr_len():
    fa = Fasta(TEST_FASTA)
    chrs, lens = fa.get_chr_len()
    assert 'chr1' in chrs
    assert lens[chrs.index('chr1')] > 0

def test_read_gc_at():
    fa = Fasta(TEST_FASTA)
    gc, at = fa.read_chromosome_gc("chr1")
    assert isinstance(gc, list)
    assert isinstance(at, list)
    assert len(gc) == len(at)

def test_ioerror_logging_and_exit(caplog):
    with caplog.at_level("ERROR"):
        with pytest.raises(SystemExit):
            Fasta("nonexistent_file.fa.gz")

    assert "Problem opening file" in caplog.text

def test_valueerror_logging_and_exit(caplog):
    with patch("cnvpytor.fasta.pysam.FastaFile", side_effect=ValueError("Fake missing index")):
        with caplog.at_level("ERROR"):
            with pytest.raises(SystemExit):
                Fasta("test_data/noindex_file.fa.gz")
    assert "Index for filename" in caplog.text

def test_detect_genome_mock(caplog):
    with patch("cnvpytor.fasta.Genome.detect_genome", return_value="hg19"):
        with caplog.at_level("INFO"):
            fa = Fasta("test_data/test.fa.gz")
    assert "Detected reference genome: hg19" in caplog.text

def test_read_gc_missing_chr(caplog):
    fasta = Fasta("test_data/test.fa.gz")

    with caplog.at_level("WARNING"):
        gc, at = fasta.read_chromosome_gc("fake_chr999")

    assert gc is None
    assert at is None
    assert "Can not find chromosome" in caplog.text

def test_gc_append_padding():
    fa = Fasta("test_data/test.fa.gz")
    with patch.object(fa, 'len', {"chr1": 1000}):  # 假設 chr1 長度為 1000
        gc, at = fa.read_chromosome_gc("chr1")
        assert gc[-1] == 0  # 最後一個 bin 是補上去

def test_p_mask_chr_not_found(caplog):
    fa = Fasta("test_data/test.fa.gz")
    with caplog.at_level("WARNING"):
        p = fa.read_chromosome_mask_p_regions("nonexist_chr")
    assert p is None
    assert "Can not find chromosome" in caplog.text

def test_p_mask_none_found():
    fa = Fasta("test_data/test.fa.gz")
    result = fa.read_chromosome_mask_p_regions("chr1")
    assert isinstance(result, list)
    assert len(result) == 0  # 沒有 P 就回傳空 list

def test_p_mask_with_span():
    fa = Fasta("test_data/test_with_p.fa.gz")
    spans = fa.read_chromosome_mask_p_regions("chr1")
    assert isinstance(spans, list)
    assert all(isinstance(x, tuple) and len(x) == 2 for x in spans)
    assert spans[0][1] > spans[0][0]

def test_print_template_four_chrs(capsys):
    fa = Fasta("test_data/test.fa.gz")
    
    mock_file = MagicMock()
    mock_file.references = ["chr1", "chr2", "chr3", "chr4"]
    mock_file.lengths = [100, 200, 300, 400]
    fa.file = mock_file  # 👈 直接用假的 file

    fa.print_reference_genome_template()

    captured = capsys.readouterr()
    assert "OrderedDict" in captured.out
    assert "\n" in captured.out