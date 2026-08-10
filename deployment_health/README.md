# Django deployment health endpoint

Register the centrally managed health URL module once in the application's
root `urls.py`:

```python
from django.urls import include, path

urlpatterns = [
    path("health/", include("deployment_health.urls")),
    # Application-owned routes follow.
]
```

Do not add `deployment_health` to `INSTALLED_APPS`: it has no models, templates
or application startup hooks. Compose and the deployment smoke test require
`GET /health/` to return a successful response.
