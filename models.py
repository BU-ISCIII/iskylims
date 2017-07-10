import datetime

from django.db import models
from django.utils import timezone
from django.utils.encoding import python_2_unicode_compatible

# Create your models here.


class Question(models.Model):
    question_text = models.CharField(max_length=200)
    pub_date = models.DateTimeField('date published')
    def __str__(self):
        return self.question_text
    
    def was_published_recently(self):
        return self.pub_date >= timezone.now() - datetime.timedelta(days=1)


class Choice(models.Model):
    question = models.ForeignKey(Question, on_delete=models.CASCADE)
    choice_text = models.CharField(max_length=200)
    votes = models.IntegerField(default=0)
    def __str__(self):
        return self.choice_text

class luis(models.Model):
	votes = models.IntegerField(default=0)
	choice_text = models.CharField(max_length=200)

########## class for upload field
class Document(models.Model):
    description = models.CharField(max_length=255, blank=True)
    csv_file = models.FileField(upload_to='documents/')
    uploaded_at = models.DateTimeField(auto_now_add=True)
    name = models.CharField(max_length=255, blank=True)
    user_id = models.CharField(max_length=50,blank=True)
    email = models.EmailField(max_length=255,blank=True)
    convert = models.BooleanField(default=False)

	
