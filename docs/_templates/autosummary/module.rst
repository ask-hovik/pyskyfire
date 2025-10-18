{{ fullname }}
{{ underline }}

.. automodule:: {{ fullname }}
   :members:
   :undoc-members:
   :inherited-members:
   :show-inheritance:

{% if modules %}
Submodules
----------

.. autosummary::
   :toctree:
   :recursive:
{% for item in modules %}
   {{ item }}
{% endfor %}
{% endif %}
