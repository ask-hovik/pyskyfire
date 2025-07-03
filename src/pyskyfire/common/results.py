import cloudpickle as pickle
from collections.abc import MutableMapping

class Results(MutableMapping):
    """A mapping you can also dot-access, with save/load built in."""
    def __init__(self):
        self._store = {}
        self.metadata = {}

    # Mapping interface
    def __getitem__(self, key):
        return self._store[key]
    def __setitem__(self, key, value):
        self._store[key] = value
    def __delitem__(self, key):
        del self._store[key]
    def __iter__(self):
        return iter(self._store)
    def __len__(self):
        return len(self._store)

    # dot-access
    def __getattr__(self, name):
        # grab the internal store *without* triggering __getattr__ again
        try:
            store = object.__getattribute__(self, "_store")
        except AttributeError:
            # _store really isn't there yet → no custom attribute
            raise AttributeError(f"{name!r} not found")
        if name in store:
            return store[name]
        raise AttributeError(f"{name!r} not found in Results")

    def __setattr__(self, name, value):
        if name in ("_store", "metadata"):
            # internal attributes go straight into __dict__
            object.__setattr__(self, name, value)
        else:
            # everyone else goes into the store
            object.__getattribute__(self, "_store")[name] = value

    # save/load
    def save(self, path):
        with open(path, "wb") as f:
            pickle.dump(self, f)
    @classmethod
    def load(cls, path):
        with open(path, "rb") as f:
            return pickle.load(f)

    # convenience
    def add(self, name, obj):
        """Explicitly register something new."""
        self._store[name] = obj
    def list_items(self):
        return list(self._store.keys())