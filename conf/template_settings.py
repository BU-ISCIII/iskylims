"""Django settings template for iSkyLIMS.

Keep application behavior and the exact Django application list in this file.
The BU-ISCIII deployment renderer replaces environment-specific database,
email, host, CSRF and secret values from the selected installation settings.
"""

import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# The renderer replaces this complete line and preserves the generated secret
# during upgrades. Never put a real production secret in this repository.
SECRET_KEY = "PLACEHOLDER"
DEBUG = djangodebug
ALLOWED_HOSTS = [
    host.strip() for host in "djangoallowedhosts".split(",") if host.strip()
]
CSRF_TRUSTED_ORIGINS = [
    origin.strip() for origin in "djangocsrftrustedorigins".split(",") if origin.strip()
]

# iSkyLIMS local applications. Add new applications here when their models,
# URLs, templates, signals or management commands must be registered by Django.
# Prefer "package.apps.AppConfigClass" when an application defines AppConfig.
INSTALLED_APPS = [
    "core",
    # "clinic",  # Enable only when the clinic application is deployed.
    "wetlab",
    "drylab",
    "django_utils",
    # Third-party applications required by iSkyLIMS.
    "mptt",
    "crispy_forms",
    "crispy_bootstrap5",
    "django_crontab",
    "django_mptt_admin",
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    "django_extensions",
    "rest_framework",
    "drf_yasg",
    "django_cleanup",
]

# Application names shown by the iSkyLIMS interface. When adding an application
# to INSTALLED_APPS, add it here only if it must appear in that interface.
APPS_NAMES = [
    ["wetlab", "Genomics unit: massive sequencing"],
    ["drylab", "Bioinformatics unit: analysis requests"],
]

MIDDLEWARE = [
    "django.middleware.security.SecurityMiddleware",
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

ROOT_URLCONF = "iskylims.urls"
WSGI_APPLICATION = "iskylims.wsgi.application"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        # iSkyLIMS stores user-managed dry-lab service templates persistently.
        "DIRS": [BASE_DIR + "/documents/drylab/services_templates"],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
                "django.template.context_processors.i18n",
            ],
        },
    }
]

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.mysql",
        "USER": "djangouser",
        "PASSWORD": "djangopass",
        "PORT": "djangoport",
        "NAME": "djangodbname",
        "HOST": os.getenv("DB_HOST", "djangohost"),
        "CONN_MAX_AGE": dbconnmaxage,
        "TEST": {
            "NAME": "iSkyLIMS_test",
        },
    }
}

AUTH_PASSWORD_VALIDATORS = [
    {
        "NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.MinimumLengthValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.CommonPasswordValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.NumericPasswordValidator",
    },
]

# Swagger currently uses HTTP Basic authentication. To use authorization tokens,
# extend SECURITY_DEFINITIONS with an apiKey entry and test schema access.
SWAGGER_SETTINGS = {"SECURITY_DEFINITIONS": {"basic": {"type": "basic"}}}

LANGUAGE_CODE = "en-us"
TIME_ZONE = "Europe/Madrid"
USE_I18N = True
USE_L10N = True
USE_TZ = False

STATIC_URL = "/static/"
STATIC_ROOT = os.path.join(BASE_DIR, "static/")
MEDIA_URL = "/documents/"
MEDIA_ROOT = os.path.join(BASE_DIR, "documents/")

CRISPY_ALLOWED_TEMPLATE_PACKS = "bootstrap5"
CRISPY_TEMPLATE_PACK = "bootstrap5"
LOGIN_REDIRECT_URL = "/"

EMAIL_BACKEND = "django.core.mail.backends.smtp.EmailBackend"
EMAIL_HOST = "emailhostserver"
EMAIL_PORT = "emailport"
EMAIL_HOST_USER = "emailhostuser"
EMAIL_HOST_PASSWORD = "emailhostpassword"
EMAIL_USE_TLS = emailhosttls
ALLOWED_EMAIL_DOMAINS = ["isciii.es", "externos.isciii.es"]

# django-crontab writes into the persistent application log directory. Add new
# jobs as (schedule, callable, output redirection) tuples and verify them with
# `python manage.py crontab show` after deployment.
LOG_CRONTAB_FILE = os.path.join(BASE_DIR, "logs", "crontab.log")
LOG_CLEAN_FILE = os.path.join(BASE_DIR, "logs", "crontab_cleanup.log")
CRONJOBS = [
    ("*/15 * * * *", "wetlab.cron.looking_for_new_runs", ">>" + LOG_CRONTAB_FILE),
]
CRONTAB_COMMAND_SUFFIX = "2>&1"

# Maximum request body retained in memory before Django streams to disk.
DATA_UPLOAD_MAX_MEMORY_SIZE = 10_000_000

# Trust this header only because Apache overwrites X-Forwarded-Proto before
# forwarding requests to Django. Never expose Gunicorn directly to untrusted clients.
SECURE_PROXY_SSL_HEADER = ("HTTP_X_FORWARDED_PROTO", "https")

DEFAULT_AUTO_FIELD = "django.db.models.AutoField"
