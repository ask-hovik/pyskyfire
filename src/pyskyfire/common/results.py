import cloudpickle as pickle
from collections.abc import MutableMapping

class Results(MutableMapping):
    """Dictionary-like container with dot access and persistence helpers.

    Behaves as a normal mutable mapping while allowing attribute-style
    access (``res.foo`` ≡ ``res["foo"]``) and built-in save/load using
    :mod:`cloudpickle`.

    Attributes
    ----------
    _store : dict[str, Any]
        Internal data store for all key/value pairs.
    metadata : dict[str, Any]
        Optional auxiliary metadata (e.g., creation date, simulation info).

    See Also
    --------
    dict : Built-in Python mapping type.
    cloudpickle : Used for serialization in :meth:`save` and :meth:`load`.

    Notes
    -----
    Dot assignment automatically writes into :attr:`_store`. Internal
    attributes such as ``_store`` and ``metadata`` are protected from
    accidental overwrite.
    """
    def __init__(self):
        self._store = {}
        self.metadata = {}

    # Mapping interface
    def __getitem__(self, key):
        """Return the value associated with ``key``."""
        return self._store[key]
    def __setitem__(self, key, value):
        """Store ``value`` under ``key``."""
        self._store[key] = value
    def __delitem__(self, key):
        """Delete ``key`` from the store."""
        del self._store[key]
    def __iter__(self):
        """Iterate over stored keys."""
        return iter(self._store)
    def __len__(self):
        """Return the number of stored items."""
        return len(self._store)

    def __getattr__(self, name):
        """Return an attribute or stored item.

        Parameters
        ----------
        name : str
            Attribute or key name.

        Returns
        -------
        Any
            The stored object corresponding to ``name``.

        Raises
        ------
        AttributeError
            If ``name`` is not found in :attr:`_store` or the instance.
        """
        try:
            store = object.__getattribute__(self, "_store")
        except AttributeError:
            # _store really isn't there yet → no custom attribute
            raise AttributeError(f"{name!r} not found")
        if name in store:
            return store[name]
        raise AttributeError(f"{name!r} not found in Results")

    def __setattr__(self, name, value):
        """Assign attributes to the internal store by default.

        Parameters
        ----------
        name : str
            Attribute name.
        value : Any
            Value to store or assign.

        Notes
        -----
        ``_store`` and ``metadata`` are written directly to ``__dict__``;
        all other names go into :attr:`_store`.
        """
        if name in ("_store", "metadata"):
            # internal attributes go straight into __dict__
            object.__setattr__(self, name, value)
        else:
            # everyone else goes into the store
            object.__getattribute__(self, "_store")[name] = value

    # save/load
    def save(self, path):
        """Serialize this :class:`Results` object to disk.

        Parameters
        ----------
        path : str | PathLike
            File path to write the pickled object.

        Notes
        -----
        Uses :func:`cloudpickle.dump` for robust serialization of arbitrary
        Python objects.
        """
        with open(path, "wb") as f:
            pickle.dump(self, f)

    @classmethod
    def load(cls, path):
        """Load a previously saved :class:`Results` object.

        Parameters
        ----------
        path : str | PathLike
            Path to the pickled file.

        Returns
        -------
        Results
            The deserialized object instance.

        Notes
        -----
        Uses :func:`cloudpickle.load` for unpickling.
        """
        with open(path, "rb") as f:
            return pickle.load(f)

    # convenience
    def add(self, name, obj):
        """Explicitly register a new object.

        Parameters
        ----------
        name : str
            Key under which to store the object.
        obj : Any
            The object to add.
        """
        self._store[name] = obj

    def list_items(self):
        """List all stored keys.

        Returns
        -------
        list[str]
            Keys currently stored in the container.
        """
        return list(self._store.keys())