"""Compatibility shims for the installed openms-insight release.

openms-insight 0.4.0 calls ``_get_render_key()`` (and, for actions sent by a
view, ``_handle_event()``) on every rendered component, but ``SequenceView``
does not subclass ``BaseComponent`` and lacks both, so rendering it raises
``AttributeError``. Imported by ``app.py`` before any page runs. Remove once a
release with t0mdavid-m/OpenMS-Insight#25 is required.
"""
from openms_insight import SequenceView

if not hasattr(SequenceView, "_get_render_key"):
    SequenceView._get_render_key = lambda self: None
if not hasattr(SequenceView, "_handle_event"):
    SequenceView._handle_event = lambda self, event: False
