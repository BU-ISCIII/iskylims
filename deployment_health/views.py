"""Minimal health views used by container orchestration and smoke tests."""

from django.http import JsonResponse


def health_check(_request):
    """Report HTTP process readiness without querying external dependencies."""
    return JsonResponse({"status": "ok"})
