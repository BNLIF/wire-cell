## Usage

```bash
python dump_json.py [filename] [alg1 alg2 ...]
```

The currently available Wire-Cell algorithms are `simple`, `charge`, `true`, `deblob`, `deadarea`, `flash` and `cluster`.
In the end, a `to_upload.zip` file is created and can be drag-to-upload to [BEE](http://www.phy.bnl.gov/wire-cell/bee/)

Some notes:

- The default algorithms to run are `simple`, `charge`, `truth` and `mc` if not specified
- The script expects the ROOT file ending with `_{eventID}.root`,
and will actually run through all files in the same directory that have the same prefix.
- If your files are in the same directory but have different prefix,
you can use a special syntax `python dump_json.py 'yourdir/*_*.root' [alg1 alg2 ...]`. Please don't foget the single quotes.

