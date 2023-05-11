"""
Django settings for iSkyLIMS project.

Generated by 'django-admin startproject' using Django 1.11.4.

For more information on this file, see
https://docs.djangoproject.com/en/1.11/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/1.11/ref/settings/
"""

import os

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/1.11/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = SECRET

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

ALLOWED_HOSTS = ["localhost", "127.0.0.1", "localserverip"]

# Application definition

INSTALLED_APPS = [
    "core",
    # "clinic",
    "wetlab",
    "drylab",
    "django_utils",
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

APPS_NAMES = [
    # ["clinic", "Clinic"],
    ["wetlab", "Masive Sequencing"],
    ["drylab", "Request Service"],
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

ROOT_URLCONF = "iSkyLIMS.urls"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
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
    },
]

WSGI_APPLICATION = "iSkyLIMS.wsgi.application"


# Database
# https://docs.djangoproject.com/en/1.11/ref/settings/#databases
DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.mysql",
        "USER": "djangouser",
        "PASSWORD": "djangopass",
        "PORT": "djangoport",
        "NAME": "djangodbname",
        "HOST": "djangohost",
        "TEST": {
            "NAME": "iSkyLIMS_test",
        },
    },
}


# Password validation
# https://docs.djangoproject.com/en/1.11/ref/settings/#auth-password-validators

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
SWAGGER_SETTINGS = {"SECURITY_DEFINITIONS": {"basic": {"type": "basic"}}}

""" For using token in the authorization request
    'api_key': {
            'type': 'apiKey',
            'in': 'header',
            'name': 'Authorization'
        }   
"""
# 'PERSIST_AUTH': True

# Internationalization
# https://docs.djangoproject.com/en/1.11/topics/i18n/

LANGUAGE_CODE = "en-us"
TIME_ZONE = "Europe/Madrid"

USE_I18N = True

USE_L10N = True

USE_TZ = False


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/1.11/howto/static-files/

STATIC_URL = "/static/"
STATIC_ROOT = os.path.join(BASE_DIR, "static/")

#  Media settings
MEDIA_URL = "/documents/"
MEDIA_ROOT = os.path.join(BASE_DIR, "documents/")

#  Crispy forms settings
CRISPY_ALLOWED_TEMPLATE_PACKS = "bootstrap5"
CRISPY_TEMPLATE_PACK = "bootstrap5"

# Redirect to home URL after login (Default redirects to /accounts/profile/)
LOGIN_REDIRECT_URL = "/"

EMAIL_BACKEND = (
    "django.core.mail.backends.console.EmailBackend"  # During development only
)
# EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'

# EMAIL settings
EMAIL_HOST = "emailhostserver"
EMAIL_PORT = "emailport"
EMAIL_HOST_USER = "emailhostuser"
EMAIL_HOST_PASSWORD = "emailhostpassword"
EMAIL_USE_TLS = emailhosttls
ALLOWED_EMAIL_DOMAINS = []

LOG_CRONTAB_FILE = os.path.join(BASE_DIR, "logs", "crontab.log")
LOG_CLEAN_FILE = os.path.join(BASE_DIR, "logs", "crontab_cleanup.log")


# Crontab settings
CRONJOBS = [
    ("*/15 * * * *", "wetlab.cron.looking_for_new_runs", ">>" + LOG_CRONTAB_FILE),
    ("0 0 1 * *", "wetlab.cron.delete_invalid_run", ">>" + LOG_CLEAN_FILE),
]

CRONTAB_COMMAND_SUFFIX = "2>&1"

DEFAULT_AUTO_FIELD = "django.db.models.AutoField"
DATA_UPLOAD_MAX_MEMORY_SIZE = 5000000
